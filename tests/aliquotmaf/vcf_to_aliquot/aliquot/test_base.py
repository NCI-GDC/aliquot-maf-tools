#!/usr/bin/env python3

import unittest
from functools import partial
from types import SimpleNamespace
from unittest import mock

from maflib.header import MafHeader
from maflib.sorter import MafSorter
from maflib.writer import MafWriter

from aliquotmaf.vcf_to_aliquot.aliquot import base as MOD


class ThisTestCase(unittest.TestCase):
    CLASS_OBJ = MOD.Aliquot

    def setUp(self) -> None:
        super().setUp()

        self.mocks = SimpleNamespace(
            HEADER=mock.MagicMock(spec_set=MafHeader),
            SORTER=mock.MagicMock(spec_set=MafSorter),
            WRITER=mock.MagicMock(spec_set=MafWriter),
        )

        self.mock_header = mock.MagicMock(spec_set=MafHeader)
        self.mocks.HEADER.from_defaults.return_value = self.mock_header

        self.mock_sorter = mock.MagicMock(spec_set=MafSorter)
        self.mocks.SORTER.return_value = self.mock_sorter

        self.mock_writer = mock.MagicMock(spec_set=MafWriter)
        self.mocks.WRITER.from_path.return_value = self.mock_writer

        self.input_vcf = "/path/to/input.vcf"
        self.output_maf = "/path/to/output.maf"
        self.tumor_sample_id = "tumor-sample-id"
        self.tumor_aliquot_uuid = "tumor-aliquot-uuid"
        self.tumor_submitter_id = "tumor-submitter-id"
        self.reference_fasta_index = "reference-fasta-index"
        self.maf_center = "maf-center"
        self.sequencer = "sequencer"
        self.src_vcf_uuid = "src_vcf_uuid"
        self.case_uuid = "case-uuid"
        self.normal_sample_id = "normal-sample-id"
        self.normal_aliquot_uuid = "normal-aliquot-uuid"
        self.normal_submitter_id = "normal-submitter-id"

    def _init_obj(self, is_tumor_only: bool = False) -> MOD.Aliquot:
        _class = partial(
            self.CLASS_OBJ,
            self.input_vcf,
            self.output_maf,
            self.tumor_sample_id,
            self.tumor_aliquot_uuid,
            self.tumor_submitter_id,
            self.reference_fasta_index,
            self.maf_center,
            self.sequencer,
            self.src_vcf_uuid,
            self.case_uuid,
            _maf_header=self.mocks.HEADER,
            _maf_sorter=self.mocks.SORTER,
            _maf_writer=self.mocks.WRITER,
        )
        if is_tumor_only:
            return _class(is_tumor_only=True)
        else:
            return _class(
                normal_sample_id=self.normal_sample_id,
                normal_aliquot_uuid=self.normal_aliquot_uuid,
                normal_submitter_id=self.normal_submitter_id,
                is_tumor_only=False,
            )

    def test_init_validation(self):
        obj = self._init_obj()
        assert obj

    def test_tumor_only_init_validation(self):
        with self.assertRaisesRegex(ValueError, "Normal ID given"):
            self.CLASS_OBJ(
                self.input_vcf,
                self.output_maf,
                self.tumor_sample_id,
                self.tumor_aliquot_uuid,
                self.tumor_submitter_id,
                self.reference_fasta_index,
                self.maf_center,
                self.sequencer,
                self.src_vcf_uuid,
                self.case_uuid,
                is_tumor_only=True,
                normal_sample_id=self.normal_sample_id,
                normal_aliquot_uuid=self.normal_aliquot_uuid,
                normal_submitter_id=self.normal_submitter_id,
                _maf_header=self.mocks.HEADER,
                _maf_sorter=self.mocks.SORTER,
                _maf_writer=self.mocks.WRITER,
            )

    @mock.patch('aliquotmaf.vcf_to_aliquot.aliquot.base.BarcodesAndCoordinate')
    def test_maf_header_initialized(self, patch):
        m = mock.Mock()
        patch.return_value = m
        obj = self._init_obj()
        header = obj.header
        self.mocks.HEADER.from_defaults.assert_called_once_with(
            version=self.CLASS_OBJ.VERSION,
            annotation=self.CLASS_OBJ.ANNOTATION,
            sort_order=m,
            fasta_index=self.reference_fasta_index,
        )

    def test_context_manager_closes_sorter_writer_on_exit(self):
        obj = self._init_obj()
        with obj as o:
            pass
        self.mock_sorter.close.assert_called_once_with()
        self.mock_writer.close.assert_called_once_with()


# __END__
