"""Main vcf2maf logic for spec gdc-1.0.0-aliquot"""
import urllib.parse
from operator import itemgetter

import pysam
from maflib.header import MafHeader, MafHeaderRecord
from maflib.sort_order import BarcodesAndCoordinate
from maflib.sorter import MafSorter
from maflib.validation import ValidationStringency
from maflib.writer import MafWriter

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


class Aliquot:
    VERSION = None
    ANNOTATION = None

    @classmethod
    def __tool_name__(cls):
        return cls.ANNOTATION

    def __init__(
        self,
        input_vcf: str,
        tumor_sample_id: str,
        normal_sample_id: str = None,
        is_tumor_only: bool = False,
    ):
        assert (normal_sample_id and not is_tumor_only) or (
            not normal_sample_id and is_tumor_only
        ), "Normal_sample_id given when tumor_only True"

        self.input_vcf = input_vcf
        self.tumor_sample_id = tumor_sample_id
        self.normal_sample_id = normal_sample_id
        self.is_tumor_only = is_tumor_only

        self._sorter = None
        self._header = None
        self._writer = None

        self._scheme = None
        self._columns = None
        self._colset = None

        # Iterate over annotators, use getattr() on self
        self.annotators = None
        self.filters = None

        # self.biotype_priority = load_json(self.options["biotype_priority_file"])
        self._biotype_priority = None
        # self.effect_priority = load_json(self.options["effect_priority_file"])
        self._effect_priority = None
        self._custom_enst = None
        """
        self.custom_enst = (
            load_enst(self.options["custom_enst"])
            if self.options["custom_enst"]
            else None
        )
        """

    @property
    def header(
        self,
        reference_fasta_index: str = None,
        tumor_only: bool = False,
        tumor_aliquot_uuid: str = None,
        normal_aliquot_uuid: str = None,
    ):
        """
        Sets up the maf header.
        """
        if not self._header:
            header = MafHeader.from_defaults(
                version=self.VERSION,
                annotation=self.ANNOTATION,
                sort_order=BarcodesAndCoordinate(),
                fasta_index=reference_fasta_index,
            )

            header_date = BaseRunner.get_header_date()
            header[header_date.key] = header_date

            if tumor_only:
                normal_aliquot = MafHeaderRecord(key="normal.aliquot", value="")
            else:
                normal_aliquot = MafHeaderRecord(
                    key="normal.aliquot", value=normal_aliquot_uuid,
                )
            header[normal_aliquot.key] = normal_aliquot

            tumor_aliquot = MafHeaderRecord(
                key="tumor.aliquot", value=tumor_aliquot_uuid
            )
            header[tumor_aliquot.key] = tumor_aliquot
            self._header = header
        return self._header

    @property
    def sorter(self, reference_fasta_index: str = None):
        if not self._sorter:
            self._sorter = MafSorter(
                max_objects_in_ram=100000,
                sort_order_name=BarcodesAndCoordinate.name(),
                scheme=self.header.scheme(),
                fasta_index=reference_fasta_index,
            )
        return self._sorter

    @property
    def writer(self, output_maf: str = None):
        if not self._writer:
            self._writer = MafWriter.from_path(
                path=output_maf,
                header=self.maf_header,
                validation_stringency=ValidationStringency.Strict,
            )
        return self._writer

    def setup_annotators(self):
        if self.annotators:
            pass
        else:
            raise ValueError("No annotators set.")
        # Add class finder for annotators

    def setup_extractors(self):
        # Add class finder for extractors
        pass

    def setup_filters(self):
        # Add class finder for filters
        pass

    @property
    def vcf(self, _pysam=pysam):
        if not self._vcf:
            vcf = _pysam.VariantFile(self.input_vcf)
            self._vcf = vcf
        return self._vcf

    def __enter__(self):
        self.vcf
        self.header()
        self.sorter()
        return self

    def __exit__(self, x, y, z):
        self.vcf.close()
        self.sorter.close()
        self.writer.close()
        for anno in self.annotators:
            if self.annotators[anno]:
                self.annotators[anno].shutdown()


class GDC_1_0_0_Aliquot(Aliquot):

    VERSION = "gdc-1.0.0"
    ANNOTATION = "gdc-1.0.0-aliquot"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Iterate over annotators, use getattr() on self
        self.annotators = (
            "dbsnp_priority_db",
            "reference_context",
            "cosmic_id",
            "mutation_status",
            "non_tcga_exac",
            "hotspots",
        )

        self.filters = (
            "common_in_exac",
            "gdc_blacklist",
            "normal_depth",
            "gdc_pon",
            "multiallelic",
            "nonexonic",
            "offtarget",
        )

        self._vcf = None

        # Annotators

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

    @property
    def vcf(self):
        if not self._vcf:
            vcf = pysam.VariantFile(self.input_vcf)

    def __enter__(self):
        # Open input and output files
        pass

    def __exit__(self, x, y, z):
        # Close file handlers
        self.vcf.close()
        sorter.close()
        if self.maf_writer:
            self.maf_writer.close()
        for anno in self.annotators:
            if self.annotators[anno]:
                self.annotators[anno].shutdown()
        pass

    def do_work(self):
        """Main wrapper function for running vcf2maf"""
        self.logger.info(
            "Processing input vcf {0}...".format(self.options["input_vcf"])
        )

        # Initialize the maf file
        self.setup_maf_header()

        sorter = self.setup_maf_sorter()

        # make these properties
        self._scheme = self.maf_header.scheme()
        self._columns = get_columns_from_header(self.maf_header)
        self._colset = set(self._columns)

        # Initialize vcf reader
        # Subclass this, add context management and iterator
        vcf_object = pysam.VariantFile(self.options["input_vcf"])

        # Make these attributes, set on init
        tumor_sample_id = self.options["tumor_vcf_id"]
        normal_sample_id = self.options["normal_vcf_id"]
        is_tumor_only = self.options["tumor_only"]

        try:
            # Validate samples
            # Move to new class
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

                # Extract data
                # Make this a NamedTuple
                data = self.extract(
                    tumor_sample_id,
                    normal_sample_id,
                    tumor_idx,
                    normal_idx,
                    ann_cols_format,
                    vep_key,
                    vcf_record,
                    is_tumor_only,
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

                # Add to sorter
                sorter += maf_record

            # Write
            self.logger.info("Writing {0} sorted records...".format(line))

            for record in sorter:

                self.maf_writer += record

            self.logger.info(
                "Finished writing {0} records".format(len(sorter))
            )  # FIXME: make sure this has len

        finally:
            vcf_object.close()
            sorter.close()
            if self.maf_writer:
                self.maf_writer.close()
            for anno in self.annotators:
                if self.annotators[anno]:
                    self.annotators[anno].shutdown()

        self.logger.info("Finished")

    # Annotators as properties?

    @property
    def mutation_status(self):
        if not self._mutation_status:
            self._mutation_status = Annotators.MutationStatus.setup(
                self._scheme, self.options["caller_id"]
            )
        return self._mutation_status

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

        if self.options["dbsnp_priority_db"]:
            self.annotators["dbsnp_priority_db"] = Annotators.DbSnpValidation.setup(
                self._scheme, self.options["dbsnp_priority_db"]
            )

        if self.options["cosmic_vcf"]:
            self.annotators["cosmic_id"] = Annotators.CosmicID.setup(
                self._scheme, self.options["cosmic_vcf"]
            )

        if self.options["non_tcga_exac_vcf"]:
            self.annotators["non_tcga_exac"] = Annotators.NonTcgaExac.setup(
                self._scheme, self.options["non_tcga_exac_vcf"]
            )

        if self.options["hotspot_tsv"]:
            self.annotators["hotspots"] = Annotators.Hotspot.setup(
                self._scheme, self.options["hotspot_tsv"]
            )

    def setup_filters(self):
        """
        Sets up all filter classes.
        """
        self.filters["common_in_exac"] = Filters.ExAC.setup(
            self.options["exac_freq_cutoff"]
        )

        self.filters["multiallelic"] = Filters.Multiallelic.setup()

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
