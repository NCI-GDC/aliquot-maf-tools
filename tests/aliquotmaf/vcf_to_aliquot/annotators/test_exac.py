"""
Tests for the ``aliquotmaf.annotators.NonTcgaExac`` class.
"""
from collections import OrderedDict

import pytest
from maflib.column_types import NullableFloatColumn
from maflib.record import MafColumnRecord

from aliquotmaf.annotators import NonTcgaExac
from aliquotmaf.converters.builder import get_builder

popkeys = ["AFR", "AMR", "EAS", "FIN", "NFE", "OTH", "SAS"]


@pytest.fixture
def setup_annotator():
    created = []

    def _make_annotator(scheme, source):
        curr = NonTcgaExac.setup(scheme, source)
        created.append(curr)
        return curr

    yield _make_annotator

    for record in created:
        record.shutdown()


@pytest.fixture
def test_scheme(get_test_scheme):
    vals = [
        ("nontcga_ExAC_AF", NullableFloatColumn),
        ("nontcga_ExAC_AF_Adj", NullableFloatColumn),
    ]
    for p in popkeys:
        vals.append(("nontcga_ExAC_AF_{0}".format(p), NullableFloatColumn))

    coldict = OrderedDict(vals)
    return get_test_scheme(coldict)


def test_setup_exac(test_scheme, setup_annotator, get_test_file):
    vcf_path = get_test_file("fake_exac.vcf.gz")
    annotator = setup_annotator(test_scheme, source=vcf_path)
    assert isinstance(annotator, NonTcgaExac)


def test_exac_annotator_1(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_test_vcf_record,
    get_empty_maf_record,
):
    """
    Simple SNP overlap test
    """
    vcf_path = get_test_file("fake_exac.vcf.gz")
    annotator = setup_annotator(test_scheme, source=vcf_path)

    vcf_record = get_test_vcf_record(
        chrom="chr1", pos=10, stop=10, ref="C", alleles=("C", "G"), alts=("G",)
    )
    expected = {
        "nontcga_ExAC_AF": 40 / float(70),
        "nontcga_ExAC_AF_Adj": 10 / float(200),
        "nontcga_ExAC_AF_AFR": 10 / float(10),
        "nontcga_ExAC_AF_AMR": 5 / float(10),
        "nontcga_ExAC_AF_EAS": 5 / float(10),
        "nontcga_ExAC_AF_FIN": 5 / float(10),
        "nontcga_ExAC_AF_NFE": 5 / float(10),
        "nontcga_ExAC_AF_OTH": 5 / float(10),
        "nontcga_ExAC_AF_SAS": 5 / float(10),
    }
    maf_record = annotator.annotate(get_empty_maf_record, vcf_record)

    for k in expected:
        assert maf_record[k].value == expected[k]


def test_exac_annotator_2(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_test_vcf_record,
    get_empty_maf_record,
):
    """
    Deletion overlap test
    """
    vcf_path = get_test_file("fake_exac.vcf.gz")
    annotator = setup_annotator(test_scheme, source=vcf_path)

    vcf_record = get_test_vcf_record(
        chrom="chr2", pos=10, stop=13, ref="ACTT", alleles=("ACTT", "A"), alts=("A",)
    )
    expected = {
        "nontcga_ExAC_AF": 40 / float(70),
        "nontcga_ExAC_AF_Adj": 10 / float(200),
        "nontcga_ExAC_AF_AFR": 10 / float(10),
        "nontcga_ExAC_AF_AMR": 5 / float(10),
        "nontcga_ExAC_AF_EAS": 5 / float(10),
        "nontcga_ExAC_AF_FIN": 5 / float(10),
        "nontcga_ExAC_AF_NFE": 5 / float(10),
        "nontcga_ExAC_AF_OTH": 5 / float(10),
        "nontcga_ExAC_AF_SAS": 5 / float(10),
    }
    maf_record = annotator.annotate(get_empty_maf_record, vcf_record)
    for k in expected:
        assert maf_record[k].value == expected[k]


def test_exac_annotator_3(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_test_vcf_record,
    get_empty_maf_record,
):
    """
    Deletion overlap test, different alt
    """
    vcf_path = get_test_file("fake_exac.vcf.gz")
    annotator = setup_annotator(test_scheme, source=vcf_path)

    vcf_record = get_test_vcf_record(
        chrom="chr2", pos=10, stop=13, ref="ACTT", alleles=("ACTT", "AG"), alts=("AG",)
    )
    expected = {
        "nontcga_ExAC_AF": 1 / float(70),
        "nontcga_ExAC_AF_Adj": 1 / float(200),
        "nontcga_ExAC_AF_AFR": 1 / float(10),
        "nontcga_ExAC_AF_AMR": 0 / float(10),
        "nontcga_ExAC_AF_EAS": 0 / float(10),
        "nontcga_ExAC_AF_FIN": 0 / float(10),
        "nontcga_ExAC_AF_NFE": 0 / float(10),
        "nontcga_ExAC_AF_OTH": 0 / float(10),
        "nontcga_ExAC_AF_SAS": 0 / float(10),
    }
    maf_record = annotator.annotate(get_empty_maf_record, vcf_record)
    for k in expected:
        assert maf_record[k].value == expected[k]


def test_exac_annotator_4(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_test_vcf_record,
    get_empty_maf_record,
):
    """
    SNP should not match due to alleles 
    """
    vcf_path = get_test_file("fake_exac.vcf.gz")
    annotator = setup_annotator(test_scheme, source=vcf_path)

    vcf_record = get_test_vcf_record(
        chrom="chr1", pos=10, stop=10, ref="C", alleles=("C", "T"), alts=("T",)
    )
    expected = {
        "nontcga_ExAC_AF": None,
        "nontcga_ExAC_AF_Adj": None,
        "nontcga_ExAC_AF_AFR": None,
        "nontcga_ExAC_AF_AMR": None,
        "nontcga_ExAC_AF_EAS": None,
        "nontcga_ExAC_AF_FIN": None,
        "nontcga_ExAC_AF_NFE": None,
        "nontcga_ExAC_AF_OTH": None,
        "nontcga_ExAC_AF_SAS": None,
    }
    maf_record = annotator.annotate(get_empty_maf_record, vcf_record)

    for k in expected:
        assert maf_record[k].value == expected[k]
