#!/usr/bin/env python3

import unittest
from types import SimpleNamespace
from unittest import mock

from aliquotmaf.vcf_to_aliquot import vcf as MOD


class ThisTestCase(unittest.TestCase):
    def setUp(self):
        pass


class TestVcfFile(ThisTestCase):
    def setUp(self):
        super().setUp()

        self.mocks = SimpleNamespace(pysam=mock.MagicMock(spec_set=MOD.pysam))

        self.variant_file = "mock.vcf"
        self.tumor_vcf_id = "tumor_id"
        self.normal_vcf_id = "normal_id"
        self.is_tumor_only = False

    def test_VcfFile_context_manager_calls_pysam_VariantFile(self):
        with MOD.VcfFile(
            variant_file=self.variant_file,
            tumor_vcf_id=self.tumor_vcf_id,
            normal_vcf_id=self.normal_vcf_id,
            is_tumor_only=self.is_tumor_only,
            di=self.mocks,
        ) as vcf:
            self.mocks.pysam.VariantFile.assert_called_once_with(self.variant_file)

    def test_VcfFile_iter_calls_fetch_on_pysam_object(self):
        mock_records = [mock.Mock()]
        mock_pysam_file = mock.MagicMock(spec_set=MOD.pysam.VariantFile)
        self.mocks.pysam.VariantFile.return_value = mock_pysam_file
        mock_pysam_file.fetch.return_value = mock_records
        with MOD.VcfFile(
            variant_file=self.variant_file,
            tumor_vcf_id=self.tumor_vcf_id,
            normal_vcf_id=self.normal_vcf_id,
            is_tumor_only=self.is_tumor_only,
            di=self.mocks,
        ) as vcf:
            vcf._tumor_idx = 1
            vcf._normal_idx = 2
            vcf._annotation_columns = ['foo']
            found_records = [record for record in vcf]
            vcf.pysam_file.fetch.assert_called_once_with()
        self.assertEqual(found_records[0].record, mock_records[0])
        self.mocks.pysam.VariantFile.assert_called_once_with(self.variant_file)


class TestVcfRecord(ThisTestCase):
    def setUp(self):
        super().setUp()

    def test_init_sets_attributes(self):
        attrs = {
            "record": 'mock_record',
            "index": 1,
            "normal_index": 2,
            "tumor_index": 3,
            "annotation_columns": ['foo', 'bar'],
        }

        obj = MOD.VcfRecord(**attrs)
        for k, v in attrs.items():
            with self.subTest(k=k):
                self.assertEqual(getattr(obj, k), v)


# __END__
