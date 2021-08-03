"""
Tests for the ``aliquotmaf.annotators.GnomAD_VCF`` class.
"""
import math
import os
from collections import OrderedDict, namedtuple

import pytest
from maflib.column_types import NullableFloatColumn, SequenceOfStrings
from maflib.record import MafColumnRecord

from aliquotmaf.annotators import GnomAD_VCF
from aliquotmaf.converters.builder import get_builder


@pytest.fixture
def setup_annotator():
    created = []

    def _make_annotator(scheme, source):
        curr = GnomAD_VCF.setup(scheme, source)
        created.append(curr)
        return curr

    yield _make_annotator

    for record in created:
        record.shutdown()


@pytest.fixture
def test_scheme(get_test_scheme):
    coldict = OrderedDict(
        [
            ("gnomAD_non_cancer_EAS_AF", NullableFloatColumn),
            ("gnomAD_non_cancer_AFR_AF", NullableFloatColumn),
            ("gnomAD_non_cancer_AMI_AF", NullableFloatColumn),
            ("gnomAD_non_cancer_MID_AF", NullableFloatColumn),
            ("gnomAD_non_cancer_SAS_AF", NullableFloatColumn),
            ("gnomAD_non_cancer_NFE_AF", NullableFloatColumn),
            ("gnomAD_non_cancer_AF", NullableFloatColumn),
            ("gnomAD_non_cancer_AMR_AF", NullableFloatColumn),
            ("gnomAD_non_cancer_OTH_AF", NullableFloatColumn),
            ("gnomAD_non_cancer_ASJ_AF", NullableFloatColumn),
            ("gnomAD_non_cancer_FIN_AF", NullableFloatColumn),
            ("gnomAD_non_cancer_MAX_AF_adj", NullableFloatColumn),
            ("gnomAD_non_cancer_MAX_AF_POPS_adj", SequenceOfStrings),
        ]
    )
    return get_test_scheme(coldict)


def test_setup_gnomad_vcf(test_scheme, setup_annotator, get_test_file):
    vcf_path = get_test_file("fake_noncancer_gnomad.vcf.gz")
    annotator = setup_annotator(test_scheme, vcf_path)
    assert isinstance(annotator, GnomAD_VCF)


def test_variant_with_maxAF(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_test_vcf_record,
    get_empty_maf_record,
):

    # setup annotator
    vcf_path = get_test_file("fake_noncancer_gnomad.vcf.gz")
    annotator = setup_annotator(test_scheme, vcf_path)

    # setup vcf record
    vcf_record = get_test_vcf_record(
        chrom="chr1", pos=10, stop=10, ref="T", alleles=("T", "A"), alts=("A",),
    )

    # annotate maf record
    maf_record = annotator.annotate(get_empty_maf_record, vcf_record)

    assert all(
        [
            math.isclose(maf_record["gnomAD_non_cancer_EAS_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_AFR_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_AMI_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_MID_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_SAS_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_NFE_AF"].value, 0.0),
            '{0:.6f}'.format(maf_record["gnomAD_non_cancer_AF"].value)
            == '{0:.6f}'.format(1.69742e-05),
            math.isclose(maf_record["gnomAD_non_cancer_AMR_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_OTH_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_ASJ_AF"].value, 0.0),
            '{0:.6f}'.format(maf_record["gnomAD_non_cancer_FIN_AF"].value)
            == '{0:.6f}'.format(0.000276778),
            '{0:.6f}'.format(maf_record["gnomAD_non_cancer_MAX_AF_adj"].value)
            == '{0:.6f}'.format(0.000276778),
            maf_record["gnomAD_non_cancer_MAX_AF_POPS_adj"].value == ['fin'],
        ]
    )


def test_without_maxAF(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_test_vcf_record,
    get_empty_maf_record,
):

    # setup annotator
    vcf_path = get_test_file("fake_noncancer_gnomad.vcf.gz")
    annotator = setup_annotator(test_scheme, vcf_path)

    # setup vcf record
    vcf_record = get_test_vcf_record(
        chrom="chr2", pos=10, stop=12, ref="CTA", alleles=("CTA", "C"), alts=("C",)
    )

    # annotate maf record
    maf_record = annotator.annotate(get_empty_maf_record, vcf_record)

    assert all(
        [
            math.isclose(maf_record["gnomAD_non_cancer_EAS_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_AFR_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_AMI_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_MID_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_SAS_AF"].value, 0.0),
            '{0:.6f}'.format(maf_record["gnomAD_non_cancer_NFE_AF"].value)
            == '{0:.6f}'.format(0.000685871),
            '{0:.6f}'.format(maf_record["gnomAD_non_cancer_AF"].value)
            == '{0:.6f}'.format(0.000890472),
            '{0:.6f}'.format(maf_record["gnomAD_non_cancer_AMR_AF"].value)
            == '{0:.6f}'.format(0.00595238),
            math.isclose(maf_record["gnomAD_non_cancer_OTH_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_ASJ_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_FIN_AF"].value, 0.0),
            maf_record["gnomAD_non_cancer_MAX_AF_adj"].value is None,
            maf_record["gnomAD_non_cancer_MAX_AF_POPS_adj"].value == [],
        ]
    )
