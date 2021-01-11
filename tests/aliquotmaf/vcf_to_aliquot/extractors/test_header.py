"""
Tests for functions parsing the VCF headers. 
"""
import pytest

from aliquotmaf.utils.utils import (
    assert_sample_in_header,
    extract_annotation_from_header,
)


def test_sample_in_header(get_test_vcf_header):
    vcf_record = get_test_vcf_header(samples=["TUMOR", "NORMAL"])
    tumor_idx = assert_sample_in_header(vcf_record, "TUMOR")
    assert tumor_idx == 9

    normal_idx = assert_sample_in_header(vcf_record, "NORMAL")
    assert normal_idx == 10

    with pytest.raises(AssertionError):
        fail_idx = assert_sample_in_header(vcf_record, "FAKE")

    fail_ok_idx = assert_sample_in_header(vcf_record, "FAKE", can_fail=True)
    assert fail_ok_idx is None


def test_extract_annotation_from_header(get_test_vcf_header):
    hdr_meta = {
        "key": "INFO",
        "items": [
            ("ID", "CSQ"),
            ("Number", "."),
            ("Type", "String"),
            (
                "Description",
                "Consequence annotations from Ensembl VEP. Format: ItemA|ItemB",
            ),
        ],
    }

    vcf_object = get_test_vcf_header(meta=hdr_meta)
    ann_cols_format, vep_key = extract_annotation_from_header(vcf_object, vep_key="CSQ")
    assert ann_cols_format == ["ItemA", "ItemB"]

    with pytest.raises(AssertionError):
        ann_cols_format, vep_key = extract_annotation_from_header(
            vcf_object, vep_key="CSX"
        )

    hdr_meta = {
        "key": "INFO",
        "items": [
            ("ID", "CSQ"),
            ("Number", "."),
            ("Type", "String"),
            (
                "Description",
                "Consequence annotations from Ensembl VEP. Format: ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF",
            ),
        ],
    }

    vcf_object = get_test_vcf_header(meta=hdr_meta)
    ann_cols_format, vep_key = extract_annotation_from_header(vcf_object, vep_key="CSQ")
    assert ann_cols_format == [
        "ExAC_AF",
        "ExAC_AF_Adj",
        "ExAC_AF_AFR",
        "ExAC_AF_AMR",
        "ExAC_AF_EAS",
        "ExAC_AF_FIN",
        "ExAC_AF_NFE",
        "ExAC_AF_OTH",
        "ExAC_AF_SAS",
    ]

    hdr_meta = {
        "key": "INFO",
        "items": [
            ("ID", "CSQ"),
            ("Number", "."),
            ("Type", "String"),
            (
                "Description",
                "Consequence annotations from Ensembl VEP. Format: ExAC_MAF|ExAC_ERROR",
            ),
        ],
    }

    vcf_object = get_test_vcf_header(meta=hdr_meta)
    with pytest.raises(KeyError):
        ann_cols_format, vep_key = extract_annotation_from_header(
            vcf_object, vep_key="CSQ"
        )
