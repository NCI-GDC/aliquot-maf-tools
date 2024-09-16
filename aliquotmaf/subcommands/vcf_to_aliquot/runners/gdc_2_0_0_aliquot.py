"""Main vcf2maf logic for spec gdc-2.0.0-aliquot"""
import urllib.parse

import pysam
from maflib.header import MafHeader, MafHeaderRecord
from maflib.sort_order import BarcodesAndCoordinate
from maflib.sorter import MafSorter
from maflib.validation import ValidationStringency
from maflib.writer import MafWriter

from aliquotmaf.constants import variant_callers
import aliquotmaf.annotators as Annotators
import aliquotmaf.filters as Filters
import aliquotmaf.subcommands.vcf_to_aliquot.extractors as Extractors
from aliquotmaf.converters.builder import get_builder
from aliquotmaf.converters.collection import InputCollection
from aliquotmaf.converters.formatters import (
    format_all_effects,
    format_alleles,
    format_depths,
    format_vcf_columns,
)
from aliquotmaf.converters.utils import get_columns_from_header, init_empty_maf_record
from aliquotmaf.subcommands.utils import (
    assert_sample_in_header,
    extract_annotation_from_header,
    load_enst,
    load_json,
)
from aliquotmaf.subcommands.vcf_to_aliquot.runners import BaseRunner


class GDC_2_0_0_Aliquot(BaseRunner):
    def __init__(self, options=dict()):
        super(GDC_2_0_0_Aliquot, self).__init__(options)

        # Load the resource files
        self.logger.info("Loading priority files")
        self.biotype_priority = load_json(self.options["biotype_priority_file"])
        self.effect_priority = load_json(self.options["effect_priority_file"])
        self.custom_enst = (
            load_enst(self.options["custom_enst"])
            if self.options["custom_enst"]
            else None
        )

        # Schema
        self.options["version"] = "gdc-1.0.0"
        self.options["annotation"] = "gdc-2.0.0-aliquot"

        # Annotators
        self.annotators = {
            "reference_context": None,
            "cosmic_id": None,
            "mutation_status": None,
            "hotspots": None,
            "gnomad": None,
            "entrez": None,
            "entrez_gene_id": None,
            "gnomad_noncancer": None,
        }

        # Filters
        self.filters = {
            "common_in_gnomAD": None,
            "gdc_blacklist": None,
            "normal_depth": None,
            "gdc_pon": None,
            "multiallelic": None,
            "nonexonic": None,
            "offtarget": None,
        }

    @classmethod
    def __validate_options__(cls, options):
        """Validates the tumor only stuff"""
        if options.tumor_only:
            options.normal_vcf_id = None
        else:
            if options.normal_aliquot_uuid is None:
                raise ValueError("--normal_aliquot_uuid is required")
            if options.normal_submitter_id is None:
                raise ValueError("--normal_submitter_id is required")
            if options.normal_bam_uuid is None:
                raise ValueError("--normal_bam_uuid is required")

    @classmethod
    def __add_arguments__(cls, parser):
        """Add the arguments to the parser"""
        vcf = parser.add_argument_group(title="VCF options")
        vcf.add_argument(
            "--tumor_only", action="store_true", help="Is this a tumor-only VCF?"
        )
        vcf.add_argument(
            "-t",
            "--tumor_vcf_id",
            default="TUMOR",
            help="Name of the tumor sample in the VCF",
        )
        vcf.add_argument(
            "-n",
            "--normal_vcf_id",
            default="NORMAL",
            help="Name of the normal sample in the VCF",
        )
        vcf.add_argument(
            "--caller_id",
            required=True,
            help="Name of the caller used to detect mutations",
            choices=variant_callers.astuple()
        )
        vcf.add_argument(
            "--src_vcf_uuid", required=True, help="The UUID of the src VCF file"
        )

        sample = parser.add_argument_group(title="Sample Metadata")
        sample.add_argument("--case_uuid", required=True, help="Sample case UUID")
        sample.add_argument(
            "--tumor_submitter_id",
            required=True,
            help="Tumor sample aliquot submitter ID",
        )
        sample.add_argument(
            "--tumor_aliquot_uuid", required=True, help="Tumor sample aliquot UUID"
        )
        sample.add_argument(
            "--tumor_bam_uuid", required=True, help="Tumor sample bam UUID"
        )

        sample.add_argument(
            "--normal_submitter_id", help="Normal sample aliquot submitter ID"
        )
        sample.add_argument("--normal_aliquot_uuid", help="Normal sample aliquot UUID")
        sample.add_argument("--normal_bam_uuid", help="Normal sample bam UUID")

        sample.add_argument("--sequencer", action="append", help="The sequencer used")
        sample.add_argument(
            "--maf_center", action="append", required=True, help="The sequencing center"
        )

        anno = parser.add_argument_group(title="Annotation Resources")
        anno.add_argument(
            "--biotype_priority_file", required=True, help="Biotype priority JSON"
        )
        anno.add_argument(
            "--effect_priority_file", required=True, help="Effect priority JSON"
        )
        anno.add_argument(
            "--custom_enst", default=None, help="Optional custom ENST overrides"
        )
        anno.add_argument(
            "--reference_fasta", required=True, help="Reference fasta file"
        )
        anno.add_argument(
            "--reference_fasta_index", required=True, help="Reference fasta fai file"
        )
        anno.add_argument(
            "--reference_context_size",
            type=int,
            default=5,
            help="Number of BP to add both upstream and "
            + "downstream from variant for reference context",
        )
        anno.add_argument(
            "--cosmic_vcf", default=None, help="Optional COSMIC VCF for annotating"
        )
        anno.add_argument("--hotspot_tsv", default=None, help="Optional hotspot TSV")
        # Entrez gene_id Annotator
        anno.add_argument(
            "--entrez_gene_id_json",
            default=None,
            help="Optional map of ensembl transcript IDs and symbols to entrez gene ID",
        )
        # GnomAD-noncancer Annotator
        # Kyle H commented out due to scalability and performance issues
        # anno.add_argument(
        #    "--gnomad_ref_prefix",
        #    default=None,
        #    help="prefix for gnomad reference files (feather format); prefix is added to `{chr}.feather`",
        # )
        anno.add_argument(
            "--gnomad_noncancer_vcf",
            default=None,
            help="Path to the bgzipped and tabix-indexed non-cancer gnomAD allele frequency VCF.",
        )

        filt = parser.add_argument_group(title="Filtering Options")

        filt.add_argument(
            "--gnomad_af_cutoff",
            default=0.001,
            type=float,
            help="Flag variants where the allele frequency in any gnomAD population "
            + "is greater than this value as common_in_gnomAD [0.001]",
        )
        filt.add_argument(
            "--gdc_blacklist",
            type=str,
            default=None,
            help="The file containing the blacklist tags and tumor aliquot uuids to "
            + "apply them to.",
        )
        filt.add_argument(
            "--min_n_depth",
            default=7,
            type=int,
            help="Flag variants where normal depth is <= INT as ndp [7].",
        )
        filt.add_argument(
            "--gdc_pon_vcf",
            type=str,
            default=None,
            help="The tabix-indexed panel of normals VCF for applying the gdc "
            + "pon filter",
        )
        filt.add_argument(
            "--nonexonic_intervals",
            type=str,
            default=None,
            help="Flag variants outside of this tabix-indexed bed file "
            + "as NonExonic",
        )
        filt.add_argument(
            "--target_intervals",
            action="append",
            help="Flag variants outside of these tabix-indexed bed files "
            + "as off_target. Use one or more times.",
        )

    def setup_maf_header(self):
        """
        Sets up the maf header.
        """
        self.maf_header = MafHeader.from_defaults(
            version=self.options["version"],
            annotation=self.options["annotation"],
            sort_order=BarcodesAndCoordinate(),
            fasta_index=self.options["reference_fasta_index"],
        )

        header_date = BaseRunner.get_header_date()
        self.maf_header[header_date.key] = header_date

        if not self.options["tumor_only"]:
            normal_aliquot = MafHeaderRecord(
                key="normal.aliquot",
                value=self.options["normal_aliquot_uuid"]
                if not self.options["tumor_only"]
                else "",
            )
            self.maf_header[normal_aliquot.key] = normal_aliquot

        tumor_aliquot = MafHeaderRecord(
            key="tumor.aliquot", value=self.options["tumor_aliquot_uuid"]
        )
        self.maf_header[tumor_aliquot.key] = tumor_aliquot

    def do_work(self):
        """Main wrapper function for running vcf2maf"""
        self.logger.info(
            "Processing input vcf {0}...".format(self.options["input_vcf"])
        )

        # Initialize the maf file
        self.setup_maf_header()

        sorter = MafSorter(
            max_objects_in_ram=100000,
            sort_order_name=BarcodesAndCoordinate.name(),
            scheme=self.maf_header.scheme(),
            fasta_index=self.options["reference_fasta_index"],
        )

        self._scheme = self.maf_header.scheme()
        self._columns = get_columns_from_header(self.maf_header)
        self._colset = set(self._columns)

        # Initialize vcf reader
        vcf_object = pysam.VariantFile(self.options["input_vcf"])
        tumor_sample_id = self.options["tumor_vcf_id"]
        normal_sample_id = self.options["normal_vcf_id"]
        is_tumor_only = self.options["tumor_only"]

        try:
            # Validate samples
            tumor_idx = assert_sample_in_header(
                vcf_object, self.options["tumor_vcf_id"]
            )
            normal_idx = assert_sample_in_header(
                vcf_object, self.options["normal_vcf_id"], can_fail=is_tumor_only
            )

            # extract annotation from header
            ann_cols_format, vep_key = extract_annotation_from_header(
                vcf_object, vep_key="CSQ"
            )

            # Initialize annotators
            self.setup_annotators()

            # Initialize filters
            self.setup_filters()

            # Convert
            line = 0
            for vcf_record in vcf_object.fetch():

                line += 1

                if line % 1000 == 0:
                    self.logger.info("Processed {0} records...".format(line))

                # Extract data
                data = self.extract(
                    tumor_sample_id,
                    normal_sample_id,
                    tumor_idx,
                    normal_idx,
                    ann_cols_format,
                    vep_key,
                    vcf_record,
                    is_tumor_only,
                    self.options["caller_id"]
                )

                # Skip rare occasions where VEP doesn't provide IMPACT or the consequence is ?
                if (
                    not data["selected_effect"]["IMPACT"]
                    or data["selected_effect"]["One_Consequence"] == "?"
                ):
                    self.logger.warn(
                        "Skipping record with unknown impact or consequence: {0} - {1}".format(
                            data["selected_effect"]["IMPACT"],
                            data["selected_effect"]["One_Consequence"],
                        )
                    )
                    continue

                # Transform
                maf_record = self.transform(
                    vcf_record, data, is_tumor_only, line_number=line
                )

                # badkeys = set(maf_record.keys()) - set(self._columns)
                # if badkeys is not set():

                #     raise KeyError("Unexpected keys {k}: {v} found in maf_record".format(
                #         k=badkeys, v=maf_record[list(badkeys)[0]]))

                # Add to sorter
                sorter += maf_record

            # Write
            self.logger.info("Writing {0} sorted records...".format(line))
            self.maf_writer = MafWriter.from_path(
                path=self.options["output_maf"],
                header=self.maf_header,
                validation_stringency=ValidationStringency.Strict,
            )

            counter = 0
            for record in sorter:

                counter += 1

                if counter % 1000 == 0:
                    self.logger.info("Wrote {0} records...".format(counter))

                self.maf_writer += record

            self.logger.info("Finished writing {0} records".format(counter))

        finally:
            vcf_object.close()
            sorter.close()
            if self.maf_writer:
                self.maf_writer.close()
            for anno in self.annotators:
                if self.annotators[anno]:
                    self.annotators[anno].shutdown()

        self.logger.info("Finished")

    def extract(
        self,
        tumor_sample_id,
        normal_sample_id,
        tumor_idx,
        normal_idx,
        ann_cols,
        vep_key,
        record,
        is_tumor_only,
        caller_id
    ):
        """
        Extract the VCF information needed to transform into MAF.
        """
        dic = {
            "var_allele_idx": None,
            "tumor_gt": None,
            "tumor_depths": None,
            "normal_gt": None,
            "normal_depths": None,
            "location_data": None,
            "effects": None,
            "selected_effect": None,
            "variant_class": None,
        }

        # Genotypes
        var_allele_idx = Extractors.VariantAlleleIndexExtractor.extract(
            tumor_genotype=record.samples[tumor_sample_id]
        )
        tumor_gt, tumor_depths = Extractors.GenotypeAndDepthsExtractor.extract(
            var_allele_idx=var_allele_idx,
            genotype=record.samples[tumor_sample_id],
            alleles=record.alleles,
            caller_id=caller_id
        )

        if not is_tumor_only:
            normal_gt, normal_depths = Extractors.GenotypeAndDepthsExtractor.extract(
                var_allele_idx=var_allele_idx,
                genotype=record.samples[normal_sample_id],
                alleles=record.alleles,
                caller_id=caller_id
            )
        else:
            normal_gt, normal_depths = None, None

        # Locations
        location_data = Extractors.LocationDataExtractor.extract(
            ref_allele=record.ref,
            var_allele=record.alleles[var_allele_idx],
            position=record.pos,
            alleles=record.alleles,
        )

        # Handle effects
        effects = Extractors.EffectsExtractor_102.extract(
            effect_priority=self.effect_priority,
            effect_keys=ann_cols,
            effect_list=[
                urllib.parse.unquote(i).split("|") for i in record.info[vep_key]
            ],
            var_idx=var_allele_idx,
        )

        effects, selected_effect = Extractors.SelectOneEffectExtractor.extract(
            all_effects=effects,
            effect_priority=self.effect_priority,
            biotype_priority=self.biotype_priority,
            custom_enst=self.custom_enst,
        )

        selected_effect = Extractors.PopulationFrequencyExtractor.extract(
            effect=selected_effect, var_allele=location_data["var_allele"]
        )

        # Handle variant class
        variant_class = Extractors.VariantClassExtractor.extract(
            cons=selected_effect["One_Consequence"],
            var_type=location_data["var_type"],
            inframe=location_data["inframe"],
        )

        # Make return dictionary
        dic["var_allele_idx"] = var_allele_idx
        dic["tumor_gt"] = tumor_gt
        dic["tumor_depths"] = tumor_depths
        dic["normal_gt"] = normal_gt
        dic["normal_depths"] = normal_depths
        dic["location_data"] = location_data
        dic["effects"] = format_all_effects(effects)
        dic["selected_effect"] = selected_effect
        dic["variant_class"] = variant_class
        dic["vcf_columns"] = format_vcf_columns(
            vcf_record=record,
            vep_key=vep_key,
            tumor_idx=tumor_idx,
            normal_idx=normal_idx,
        )
        return dic

    def transform(self, vcf_record, data, is_tumor_only, line_number=None):
        """
        Transform into maf record.
        """

        # Generic data
        collection = InputCollection()

        collection.add(
            column="Hugo_Symbol",
            value=data["selected_effect"].get("Hugo_Symbol"),
            default="Unknown",
        )

        collection.add(column="Center", value=self.options["maf_center"])
        collection.add(column="NCBI_Build", value="GRCh38")
        collection.add(column="Chromosome", value=vcf_record.chrom)
        collection.add(column="Start_Position", value=data["location_data"]["start"])
        collection.add(column="End_Position", value=data["location_data"]["stop"])
        collection.add(column="Strand", value="+")
        collection.add(column="Variant_Classification", value=data["variant_class"])
        collection.add(column="Variant_Type", value=data["location_data"]["var_type"])
        collection.add(
            column="Reference_Allele", value=data["location_data"]["ref_allele"]
        )

        for k, v in zip(
            ["Tumor_Seq_Allele1", "Tumor_Seq_Allele2"],
            format_alleles(
                genotype=data["tumor_gt"],
                alleles=data["location_data"]["alleles"],
                defaults=[
                    data["location_data"]["ref_allele"],
                    data["location_data"]["var_allele"],
                ],
            ),
        ):
            collection.add(column=k, value=v)

        if not is_tumor_only:
            for k, v in zip(
                ["Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2"],
                format_alleles(
                    genotype=data["normal_gt"],
                    alleles=data["location_data"]["alleles"],
                    defaults=[
                        data["location_data"]["ref_allele"],
                        data["location_data"]["ref_allele"],
                    ],
                ),
            ):
                collection.add(column=k, value=v)
        else:
            for k in ["Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2"]:
                collection.add(column=k, value="")

        collection.add(column="dbSNP_RS", value=data["selected_effect"]["dbSNP_RS"])

        collection.add(
            column="Tumor_Sample_Barcode", value=self.options["tumor_submitter_id"]
        )
        collection.add(
            column="Matched_Norm_Sample_Barcode",
            value=self.options["normal_submitter_id"],
            default="",
        )
        collection.add(column="Sequencer", value=self.options["sequencer"], default="")
        collection.add(
            column="Tumor_Sample_UUID", value=self.options["tumor_aliquot_uuid"]
        )
        collection.add(
            column="Matched_Norm_Sample_UUID",
            value=self.options["normal_aliquot_uuid"],
            default="",
        )
        collection.add(column="all_effects", value=";".join(data["effects"]))

        for k, v in zip(
            ["t_depth", "t_ref_count", "t_alt_count"],
            format_depths(
                genotype=data["tumor_gt"],
                depths=data["tumor_depths"],
                var_allele_idx=data["var_allele_idx"],
                default_total_dp=0,
            ),
        ):
            collection.add(column=k, value=v)

        if not is_tumor_only:
            for k, v in zip(
                ["n_depth", "n_ref_count", "n_alt_count"],
                format_depths(
                    genotype=data["normal_gt"],
                    depths=data["normal_depths"],
                    var_allele_idx=data["var_allele_idx"],
                ),
            ):
                collection.add(column=k, value=v)
        else:
            for k in ["n_depth", "n_ref_count", "n_alt_count"]:
                collection.add(column=k, value=None)

        # Add other columns from selected_effect if they are canonical in the schema
        for k in data["selected_effect"]:
            if k in self._colset and k not in collection._colset:
                collection.add(column=k, value=data["selected_effect"][k])

        # Set other uuids
        collection.add(column="src_vcf_id", value=self.options["src_vcf_uuid"])
        collection.add(column="tumor_bam_uuid", value=self.options["tumor_bam_uuid"])
        collection.add(column="normal_bam_uuid", value=self.options["normal_bam_uuid"])
        collection.add(column="case_id", value=self.options["case_uuid"])

        # VCF columns
        collection.add(column="FILTER", value=";".join(sorted(list(vcf_record.filter))))
        collection.add(column="vcf_region", value=data["vcf_columns"]["vcf_region"])
        collection.add(column="vcf_info", value=data["vcf_columns"]["vcf_info"])
        collection.add(column="vcf_format", value=data["vcf_columns"]["vcf_format"])
        collection.add(column="vcf_tumor_gt", value=data["vcf_columns"]["vcf_tumor_gt"])
        collection.add(
            column="vcf_normal_gt", value=data["vcf_columns"].get("vcf_normal_gt")
        )

        # Set the other columns to none
        collection.add(column="Score", value="")
        collection.add(column="BAM_File", value="")
        collection.add(column="Sequencing_Phase", value="")

        # I don't think this is needed, all of these are created at initialization from the schema columns
        # dbSNP_Val_Status needs to remain as column but will always be empty
        anno_set = ("COSMIC", "CONTEXT", "Mutation_Status")
        for i in self._colset - set(collection.columns()):
            if i not in anno_set:
                collection.add(column=i, value=None)
        collection.transform(self._scheme)

        # Generate maf record
        maf_record = init_empty_maf_record(line_number=line_number)
        for i in collection:
            maf_record += i.transformed

        # # check if None introduced before here <------ DEBUG
        # foo = set(maf_record.keys()) - set(self._columns)
        # if len(foo) != 0:
        #     raise KeyError("Unexpected keys found: {}".format(foo))

        # Annotations
        maf_record["dbSNP_Val_Status"] = get_builder(
            "dbSNP_Val_Status", self._scheme, value=None
        )

        if self.annotators["cosmic_id"]:
            maf_record = self.annotators["cosmic_id"].annotate(maf_record, vcf_record)
        else:
            maf_record["COSMIC"] = get_builder("COSMIC", self._scheme, value=None)

        if self.annotators["hotspots"]:
            maf_record = self.annotators["hotspots"].annotate(maf_record)
        else:
            maf_record["hotspot"] = get_builder("hotspot", self._scheme, value=None)

        if self.annotators["entrez_gene_id"]:
            maf_record = self.annotators["entrez_gene_id"].annotate(maf_record)
        else:
            maf_record["Entrez_Gene_Id"] = get_builder(
                "entrez_gene_id", self._scheme, value=0
            )

        if self.annotators["gnomad_noncancer"]:
            maf_record = self.annotators["gnomad_noncancer"].annotate(
                maf_record, vcf_record, data["var_allele_idx"]
            )

        maf_record = self.annotators["reference_context"].annotate(
            maf_record, vcf_record
        )
        maf_record = self.annotators["mutation_status"].annotate(
            maf_record, vcf_record, self.options["tumor_vcf_id"]
        )

        # Filters
        gdc_filters = []
        for filt_key in self.filters:
            filt_obj = self.filters[filt_key]
            if filt_obj and filt_obj.filter(maf_record):
                gdc_filters.extend(filt_obj.tags)

        maf_record["GDC_FILTER"] = get_builder(
            "GDC_FILTER", self._scheme, value=";".join(sorted(gdc_filters))
        )

        return maf_record

    def setup_annotators(self):
        """
        Sets up all annotator classes.
        """
        self.annotators["mutation_status"] = Annotators.MutationStatus.setup(
            self._scheme, self.options["caller_id"]
        )

        self.annotators["reference_context"] = Annotators.ReferenceContext.setup(
            self._scheme,
            self.options["reference_fasta"],
            self.options["reference_context_size"],
        )

        if self.options["cosmic_vcf"]:
            self.annotators["cosmic_id"] = Annotators.CosmicID.setup(
                self._scheme, self.options["cosmic_vcf"]
            )

        if self.options["hotspot_tsv"]:
            self.annotators["hotspots"] = Annotators.Hotspot.setup(
                self._scheme, self.options["hotspot_tsv"]
            )

        if self.options["entrez_gene_id_json"]:
            self.annotators["entrez_gene_id"] = Annotators.Entrez.setup(
                self._scheme, self.options["entrez_gene_id_json"]
            )

        # Commented out for now due to performance
        # if self.options["gnomad_ref_prefix"]:
        #    self.annotators["gnomad_noncancer"] = Annotators.GnomAD_VCF.setup(
        #        self._scheme, self.options["gnomad_ref_prefix"]
        #    )

        if self.options["gnomad_noncancer_vcf"]:
            self.annotators["gnomad_noncancer"] = Annotators.GnomAD_VCF.setup(
                self._scheme, self.options["gnomad_noncancer_vcf"]
            )

    def setup_filters(self):
        """
        Sets up all filter classes.
        """

        self.filters["multiallelic"] = Filters.Multiallelic.setup()

        if self.options["gnomad_af_cutoff"]:
            self.filters["common_in_gnomAD"] = Filters.FilterGnomAD.setup(
                self.options["gnomad_af_cutoff"]
            )

        if self.options["gdc_blacklist"]:
            self.filters["gdc_blacklist"] = Filters.GdcBlacklist.setup(
                self.options["gdc_blacklist"]
            )

        if not self.options["tumor_only"]:
            self.filters["normal_depth"] = Filters.NormalDepth.setup(
                self.options["min_n_depth"]
            )

        if self.options["gdc_pon_vcf"]:
            self.filters["gdc_pon"] = Filters.GdcPon.setup(self.options["gdc_pon_vcf"])

        if self.options["nonexonic_intervals"]:
            self.filters["nonexonic"] = Filters.NonExonic.setup(
                self.options["nonexonic_intervals"]
            )

        if self.options["target_intervals"]:
            self.filters["off_target"] = Filters.OffTarget.setup(
                self.options["target_intervals"]
            )

    @classmethod
    def __tool_name__(cls):
        return "gdc-2.0.0-aliquot"
