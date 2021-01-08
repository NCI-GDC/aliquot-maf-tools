#!/usr/bin/env python3
"""Base Aliquot MAF class.

Children classes should represent a single MAF spec, each requiring annotators
    and filters be specified.

Example:
    >>> with aliquot(*args, **kwargs) as aliquot_object:
    >>>     aliquot_object.setup_annotators(args)
    >>>     aliquot_object.setup_filters(args)

"""

import datetime
from typing import Optional

from maflib.header import MafHeader, MafHeaderRecord
from maflib.sort_order import BarcodesAndCoordinate
from maflib.sorter import MafSorter
from maflib.validation import ValidationStringency
from maflib.writer import MafWriter

from aliquotmaf.vcf_to_aliquot.converters.utils import get_columns_from_header


class Aliquot:
    """Base Aliquot MAF class.
    Children classes should represent a single MAF spec, each requiring annotators
        and filters be specified.

    Properties:
        ANNOTATION (list): MAF annotation
        VERSION (list): MAF version
    """

    ANNOTATION = None
    VERSION = None

    ANNOTATORS = None
    FILTERS = None

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
        _maf_header=MafHeader,
        _maf_sorter=MafSorter,
        _maf_writer=MafWriter,
    ):

        self.tumor_sample_id = tumor_sample_id
        self.tumor_aliquot_uuid = tumor_aliquot_uuid
        self.normal_sample_id = normal_sample_id
        self.normal_aliquot_uuid = normal_aliquot_uuid
        self.is_tumor_only = is_tumor_only
        self.tumor_submitter_id = tumor_submitter_id
        self.maf_center = maf_center
        self.sequencer = sequencer
        self.normal_submitter_id = normal_submitter_id
        self.src_vcf_uuid = src_vcf_uuid
        self.case_uuid = case_uuid

        # Files
        self.input_vcf = input_vcf
        self.reference_fasta_index = reference_fasta_index
        self.output_maf = output_maf

        # Mocks
        self._maf_header = _maf_header
        self._maf_sorter = _maf_sorter
        self._maf_writer = _maf_writer

        # Properties
        self._header: MafHeader = None
        self._sorter: MafSorter = None
        self._writer: MafWriter = None

        # MafHeader information
        self._scheme = None
        self._columns = None
        self._colset = None

        # To be set by children classes
        self.annotators = None
        self.filters = None

        self.validate_tumor_only()

    @property
    def header(self) -> MafHeader:
        """
        Sets up the maf header.
        """
        if not self._header:
            header = self._maf_header.from_defaults(
                version=self.VERSION,
                annotation=self.ANNOTATION,
                sort_order=BarcodesAndCoordinate(),
                fasta_index=self.reference_fasta_index,
            )

            header_date = MafHeaderRecord(
                key="filedate", value=datetime.date.today().strftime("%Y%m%d")
            )
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
            self._sorter = self._maf_sorter(
                max_objects_in_ram=100000,
                sort_order_name=BarcodesAndCoordinate.name(),
                scheme=self.header.scheme(),
                fasta_index=self.reference_fasta_index,
            )
        return self._sorter

    @property
    def writer(self) -> MafWriter:
        if not self._writer:
            self._writer = self._maf_writer.from_path(
                path=self.output_maf,
                header=self.header,
                validation_stringency=ValidationStringency.Strict,
            )
        return self._writer

    @property
    def scheme(self):
        if not self._scheme:
            self._scheme = self.header.scheme()
        return self._scheme

    @property
    def colset(self):
        if not self._colset:
            self._colset = list(set(self.columns))
        return self._colset

    @property
    def columns(self):
        if not self._columns:
            self._columns = get_columns_from_header(self.header)
        return self._columns

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.sorter.close()
        self.writer.close()
        if self.annotators:
            for anno in self.annotators.values():
                anno.shutdown()

    def setup_annotators(self):
        raise NotImplementedError()

    def setup_filters(self):
        raise NotImplementedError()

    def validate_tumor_only(self):
        """Asserts no normal information given while tumor-only, and vice versa."""
        if self.is_tumor_only and self.normal_sample_id is not None:
            msg = "Normal ID given for tumor-only aliquot."
            raise ValueError(msg)
        elif not self.is_tumor_only and not self.normal_sample_id:
            msg = "Missing Normal ID"
            raise ValueError(msg)
        else:
            pass


# __END__
