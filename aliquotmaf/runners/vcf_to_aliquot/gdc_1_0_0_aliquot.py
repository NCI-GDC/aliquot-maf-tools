"""Main vcf2maf logic for spec gdc-1.0.0-aliquot"""
import urllib.parse
from operator import itemgetter
from typing import Optional

import pysam
from maflib.column_types import MutationStatus
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
from aliquotmaf.vcf_to_aliquot.annotators.cosmic import CosmicID
from aliquotmaf.vcf_to_aliquot.annotators.dbsnp_validation import DbSnpValidation
from aliquotmaf.vcf_to_aliquot.annotators.hotspot import Hotspot
from aliquotmaf.vcf_to_aliquot.annotators.nontcga_exac import NonTcgaExac
from aliquotmaf.vcf_to_aliquot.annotators.reference_context import ReferenceContext
from aliquotmaf.vcf_to_aliquot.filters import (
    ExAC,
    GdcBlacklist,
    GdcPon,
    Multiallelic,
    NonExonic,
    NormalDepth,
    OffTarget,
)


class Aliquot:
    VERSION = None
    ANNOTATION = None

    @classmethod
    def __tool_name__(cls):
        return cls.ANNOTATION

    def __init__(
        self,
        input_vcf: str,
        output_maf: str,
        tumor_sample_id: str,
        tumor_aliquot_uuid: str,
        tumor_submitter_id: str,
        reference_fasta_index: str,
        maf_center: str,
        sequencer: str,
        src_vcf_uuid: str,
        case_uuid: str,
        normal_sample_id: Optional[str] = None,
        normal_aliquot_uuid: Optional[str] = None,
        normal_submitter_id: Optional[str] = None,
        is_tumor_only: bool = False,
    ):
        """
        assert (normal_sample_id and not is_tumor_only) or (
            not normal_sample_id and is_tumor_only
        ), "Normal_sample_id given when tumor_only True"
        """
        self.validate_tumor_only()

        self.input_vcf = input_vcf
        self.tumor_sample_id = tumor_sample_id
        self.tumor_aliquot_uuid = tumor_aliquot_uuid
        self.normal_sample_id = normal_sample_id
        self.normal_aliquot_uuid = normal_aliquot_uuid
        self.is_tumor_only = is_tumor_only
        self.reference_fasta_index = reference_fasta_index
        self.output_maf = output_maf
        self.tumor_submitter_id = tumor_submitter_id
        self.maf_center = maf_center
        self.sequencer = sequencer
        self.normal_submitter_id = normal_submitter_id
        self.src_vcf_uuid = src_vcf_uuid
        self.case_uuid = case_uuid

        self._header: MafHeader = None
        self._sorter: MafSorter = None
        self._writer: MafWriter = None

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
    def header(self) -> MafHeader:
        """
        Sets up the maf header.
        """
        if not self._header:
            header = MafHeader.from_defaults(
                version=self.VERSION,
                annotation=self.ANNOTATION,
                sort_order=BarcodesAndCoordinate(),
                fasta_index=self.reference_fasta_index,
            )

            header_date = BaseRunner.get_header_date()
            header[header_date.key] = header_date

            if self.is_tumor_only:
                normal_aliquot = MafHeaderRecord(key="normal.aliquot", value="")
            else:
                normal_aliquot = MafHeaderRecord(
                    key="normal.aliquot", value=self.normal_aliquot_uuid,
                )
            header[normal_aliquot.key] = normal_aliquot

            tumor_aliquot = MafHeaderRecord(
                key="tumor.aliquot", value=self.tumor_aliquot_uuid
            )
            header[tumor_aliquot.key] = tumor_aliquot
            self._header = header
        return self._header

    @property
    def sorter(self) -> MafSorter:
        if not self._sorter:
            self._sorter = MafSorter(
                max_objects_in_ram=100000,
                sort_order_name=BarcodesAndCoordinate.name(),
                scheme=self.header.scheme(),
                fasta_index=self.reference_fasta_index,
            )
        return self._sorter

    @property
    def writer(self) -> MafWriter:
        if not self._writer:
            self._writer = MafWriter.from_path(
                path=self.output_maf,
                header=self.header,
                validation_stringency=ValidationStringency.Strict,
            )
        return self._writer

    @property
    def scheme(self):
        if not self._scheme:
            scheme = self.header.scheme()
            self._scheme = scheme
        return self._scheme

    @property
    def colset(self):
        if not self._colset:
            self._colset = set(self.columns)
        return self._colset

    @property
    def columns(self):
        if not self._columns:
            self._columns = get_columns_from_header(self.header)
        return self._columns

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

    def __enter__(self):
        return self

    def __exit__(self, x, y, z):
        self.sorter.close()
        self.writer.close()
        for anno in self.annotators.values():
            anno.shutdown()

    def validate_tumor_only(self):
        pass


class GDC_1_0_0_Aliquot(Aliquot):

    VERSION = "gdc-1.0.0"
    ANNOTATION = "gdc-1.0.0-aliquot"

    ANNOTATORS = (
        MutationStatus,
        ReferenceContext,
        DbSnpValidation,
        CosmicID,
        NonTcgaExac,
        Hotspot,
    )
    FILTERS = (
        ExAC,
        GdcBlacklist,
        NormalDepth,
        GdcPon,
        Multiallelic,
        NonExonic,
        OffTarget,
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

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
            # Moved to new VCF class
            tumor_idx = assert_sample_in_header(
                vcf_object, self.options["tumor_vcf_id"]
            )
            # Moved to new VCF class
            normal_idx = assert_sample_in_header(
                vcf_object, self.options["normal_vcf_id"], can_fail=is_tumor_only
            )

            # extract annotation from header
            # Moved to new VCF class
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

    def setup_annotators(self, args):
        """
        Sets up all annotator classes.
        """
        if not self.annotators:
            self.annotators = {}
            for Annotator in self.ANNOTATORS:
                self.annotators[Annotator.__name__] = Annotator.setup(self.scheme, args)

        """
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
        """

    def setup_filters(self, args):
        """
        Sets up all filter classes.
        """
        if not self.filters:
            self.filters = {}
            for Filter in self.FILTERS:
                filt = Filter.setup(args)
                if filt:
                    self.filters[Filter.__name__] = filt


# __END__
