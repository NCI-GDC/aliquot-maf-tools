"""
Tests for the ``aliquotmaf.annotators.GnomAD`` class.
"""
import math
from collections import OrderedDict, namedtuple

import pytest
from maflib.column_types import NullableFloatColumn, SequenceOfStrings
from maflib.record import MafColumnRecord

from aliquotmaf.annotators import GnomAD
from aliquotmaf.converters.builder import get_builder


@pytest.fixture
def setup_annotator():
    created = []

    def _make_annotator(scheme, refpath, refpattern):
        curr = GnomAD.setup(scheme, refpath, refpattern)
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


def test_setup_gnomad(test_scheme, setup_annotator, get_test_file):
    ref_path = get_test_file("gnomad")
    ref_pattern = 'gnomad_test.{}.feather'
    annotator = setup_annotator(test_scheme, ref_path, ref_pattern)
    assert isinstance(annotator, GnomAD)


def test_variant_with_maxAF(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_test_vcf_record,
    get_empty_maf_record,
):

    # setup annotator
    ref_path = get_test_file("gnomad")
    ref_pattern = 'gnomad_test.{}.feather'
    annotator = setup_annotator(test_scheme, ref_path, ref_pattern)

    # setup vcf record
    vcf_record = get_test_vcf_record(
        chrom="chr1",
        pos=10123,
        stop=10123,
        ref="CCCTAA",
        alleles=("CCCTAA", "C"),
        alts=("C",),
    )

    # annotate maf record
    maf_record = annotator.annotate(get_empty_maf_record, vcf_record)

    print(
        [
            maf_record["gnomAD_non_cancer_EAS_AF"].value,
            maf_record["gnomAD_non_cancer_AFR_AF"].value,
            maf_record["gnomAD_non_cancer_AMI_AF"].value,
            maf_record["gnomAD_non_cancer_MID_AF"].value,
            maf_record["gnomAD_non_cancer_SAS_AF"].value,
            maf_record["gnomAD_non_cancer_NFE_AF"].value,
            maf_record["gnomAD_non_cancer_AF"].value,
            maf_record["gnomAD_non_cancer_AMR_AF"].value,
            maf_record["gnomAD_non_cancer_OTH_AF"].value,
            maf_record["gnomAD_non_cancer_ASJ_AF"].value,
            maf_record["gnomAD_non_cancer_FIN_AF"].value,
            maf_record["gnomAD_non_cancer_MAX_AF_adj"].value,
            maf_record["gnomAD_non_cancer_MAX_AF_POPS_adj"].value,
        ]
    )

    assert all(
        [
            math.isclose(maf_record["gnomAD_non_cancer_EAS_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_AFR_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_AMI_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_MID_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_SAS_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_NFE_AF"].value, 5.47106e-05),
            math.isclose(maf_record["gnomAD_non_cancer_AF"].value, 2.64894e-05),
            math.isclose(maf_record["gnomAD_non_cancer_AMR_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_OTH_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_ASJ_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_FIN_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_MAX_AF_adj"].value, 5.47106e-05),
            maf_record["gnomAD_non_cancer_MAX_AF_POPS_adj"].value == ['nfe'],
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
    ref_path = get_test_file("gnomad")
    ref_pattern = 'gnomad_test.{}.feather'
    annotator = setup_annotator(test_scheme, ref_path, ref_pattern)

    # setup vcf record
    vcf_record = get_test_vcf_record(
        chrom="chr1", pos=10111, stop=10111, ref="C", alleles=("C", "A"), alts=("A",)
    )

    # annotate maf record
    maf_record = annotator.annotate(get_empty_maf_record, vcf_record)

    print(
        [
            maf_record["gnomAD_non_cancer_EAS_AF"].value,
            maf_record["gnomAD_non_cancer_AFR_AF"].value,
            maf_record["gnomAD_non_cancer_AMI_AF"].value,
            maf_record["gnomAD_non_cancer_MID_AF"].value,
            maf_record["gnomAD_non_cancer_SAS_AF"].value,
            maf_record["gnomAD_non_cancer_NFE_AF"].value,
            maf_record["gnomAD_non_cancer_AF"].value,
            maf_record["gnomAD_non_cancer_AMR_AF"].value,
            maf_record["gnomAD_non_cancer_OTH_AF"].value,
            maf_record["gnomAD_non_cancer_ASJ_AF"].value,
            maf_record["gnomAD_non_cancer_FIN_AF"].value,
            maf_record["gnomAD_non_cancer_MAX_AF_adj"].value,
            maf_record["gnomAD_non_cancer_MAX_AF_POPS_adj"].value,
        ]
    )

    assert all(
        [
            math.isclose(maf_record["gnomAD_non_cancer_EAS_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_AFR_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_AMI_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_MID_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_SAS_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_NFE_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_AF"].value, 2.28092e-05),
            math.isclose(maf_record["gnomAD_non_cancer_AMR_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_OTH_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_ASJ_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_FIN_AF"].value, 0.000486855),
            maf_record["gnomAD_non_cancer_MAX_AF_adj"].value is None,
            maf_record["gnomAD_non_cancer_MAX_AF_POPS_adj"].value == [],
        ]
    )


def test_multiple_variants_at_position(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_test_vcf_record,
    get_empty_maf_record,
):

    # setup annotator
    ref_path = get_test_file("gnomad")
    ref_pattern = 'gnomad_test.{}.feather'
    annotator = setup_annotator(test_scheme, ref_path, ref_pattern)

    # setup vcf record
    vcf_record = get_test_vcf_record(
        chrom="chr1", pos=10126, stop=10126, ref="T", alleles=("T", "G"), alts=("G",)
    )

    # annotate maf record
    maf_record = annotator.annotate(get_empty_maf_record, vcf_record)

    print(
        [
            maf_record["gnomAD_non_cancer_EAS_AF"].value,
            maf_record["gnomAD_non_cancer_AFR_AF"].value,
            maf_record["gnomAD_non_cancer_AMI_AF"].value,
            maf_record["gnomAD_non_cancer_MID_AF"].value,
            maf_record["gnomAD_non_cancer_SAS_AF"].value,
            maf_record["gnomAD_non_cancer_NFE_AF"].value,
            maf_record["gnomAD_non_cancer_AF"].value,
            maf_record["gnomAD_non_cancer_AMR_AF"].value,
            maf_record["gnomAD_non_cancer_OTH_AF"].value,
            maf_record["gnomAD_non_cancer_ASJ_AF"].value,
            maf_record["gnomAD_non_cancer_FIN_AF"].value,
            maf_record["gnomAD_non_cancer_MAX_AF_adj"].value,
            maf_record["gnomAD_non_cancer_MAX_AF_POPS_adj"].value,
        ]
    )

    assert all(
        [
            math.isclose(maf_record["gnomAD_non_cancer_EAS_AF"].value, 0.000974659),
            math.isclose(maf_record["gnomAD_non_cancer_AFR_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_AMI_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_MID_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_SAS_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_NFE_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_AF"].value, 2.37744e-05),
            math.isclose(maf_record["gnomAD_non_cancer_AMR_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_OTH_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_ASJ_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_FIN_AF"].value, 0.0),
            maf_record["gnomAD_non_cancer_MAX_AF_adj"].value is None,
            maf_record["gnomAD_non_cancer_MAX_AF_POPS_adj"].value == [],
        ]
    )


def test_variant_without_annotation(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_test_vcf_record,
    get_empty_maf_record,
):

    # setup annotator
    ref_path = get_test_file("gnomad")
    ref_pattern = 'gnomad_test.{}.feather'
    annotator = setup_annotator(test_scheme, ref_path, ref_pattern)

    # setup vcf record
    vcf_record = get_test_vcf_record(
        chrom="chr1", pos=10112, stop=10112, ref="C", alleles=("C", "A"), alts=("A",)
    )

    # annotate maf record
    maf_record = annotator.annotate(get_empty_maf_record, vcf_record)

    print(
        [
            maf_record["gnomAD_non_cancer_EAS_AF"].value,
            maf_record["gnomAD_non_cancer_AFR_AF"].value,
            maf_record["gnomAD_non_cancer_AMI_AF"].value,
            maf_record["gnomAD_non_cancer_MID_AF"].value,
            maf_record["gnomAD_non_cancer_SAS_AF"].value,
            maf_record["gnomAD_non_cancer_NFE_AF"].value,
            maf_record["gnomAD_non_cancer_AF"].value,
            maf_record["gnomAD_non_cancer_AMR_AF"].value,
            maf_record["gnomAD_non_cancer_OTH_AF"].value,
            maf_record["gnomAD_non_cancer_ASJ_AF"].value,
            maf_record["gnomAD_non_cancer_FIN_AF"].value,
            maf_record["gnomAD_non_cancer_MAX_AF_adj"].value,
            maf_record["gnomAD_non_cancer_MAX_AF_POPS_adj"].value,
        ]
    )

    assert all(
        [
            maf_record["gnomAD_non_cancer_EAS_AF"].value is None,
            maf_record["gnomAD_non_cancer_AFR_AF"].value is None,
            maf_record["gnomAD_non_cancer_AMI_AF"].value is None,
            maf_record["gnomAD_non_cancer_MID_AF"].value is None,
            maf_record["gnomAD_non_cancer_SAS_AF"].value is None,
            maf_record["gnomAD_non_cancer_NFE_AF"].value is None,
            maf_record["gnomAD_non_cancer_AF"].value is None,
            maf_record["gnomAD_non_cancer_AMR_AF"].value is None,
            maf_record["gnomAD_non_cancer_OTH_AF"].value is None,
            maf_record["gnomAD_non_cancer_ASJ_AF"].value is None,
            maf_record["gnomAD_non_cancer_FIN_AF"].value is None,
            maf_record["gnomAD_non_cancer_MAX_AF_adj"].value is None,
            maf_record["gnomAD_non_cancer_MAX_AF_POPS_adj"].value == [],
        ]
    )


def test_traverse_chromosomes(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_test_vcf_record,
    get_empty_maf_record,
):

    # setup annotator
    ref_path = get_test_file("gnomad")
    ref_pattern = 'gnomad_test.{}.feather'
    annotator = setup_annotator(test_scheme, ref_path, ref_pattern)

    # record on chr1
    vcf_record = get_test_vcf_record(
        chrom="chr1", pos=10126, stop=10126, ref="T", alleles=("T", "G"), alts=("G",)
    )
    maf_record = annotator.annotate(get_empty_maf_record, vcf_record)

    # record on chr4
    vcf_record = get_test_vcf_record(
        chrom="chr4", pos=10147, stop=10147, ref="C", alleles=("C", "G"), alts=("G",)
    )
    maf_record = annotator.annotate(get_empty_maf_record, vcf_record)

    print(
        [
            maf_record["gnomAD_non_cancer_EAS_AF"].value,
            maf_record["gnomAD_non_cancer_AFR_AF"].value,
            maf_record["gnomAD_non_cancer_AMI_AF"].value,
            maf_record["gnomAD_non_cancer_MID_AF"].value,
            maf_record["gnomAD_non_cancer_SAS_AF"].value,
            maf_record["gnomAD_non_cancer_NFE_AF"].value,
            maf_record["gnomAD_non_cancer_AF"].value,
            maf_record["gnomAD_non_cancer_AMR_AF"].value,
            maf_record["gnomAD_non_cancer_OTH_AF"].value,
            maf_record["gnomAD_non_cancer_ASJ_AF"].value,
            maf_record["gnomAD_non_cancer_FIN_AF"].value,
            maf_record["gnomAD_non_cancer_MAX_AF_adj"].value,
            maf_record["gnomAD_non_cancer_MAX_AF_POPS_adj"].value,
        ]
    )

    assert all(
        [
            math.isclose(maf_record["gnomAD_non_cancer_EAS_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_AFR_AF"].value, 7.08065e-05),
            math.isclose(maf_record["gnomAD_non_cancer_AMI_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_MID_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_SAS_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_NFE_AF"].value, 7.81779e-05),
            math.isclose(maf_record["gnomAD_non_cancer_AF"].value, 7.7738e-05),
            math.isclose(maf_record["gnomAD_non_cancer_AMR_AF"].value, 0.000229253),
            math.isclose(maf_record["gnomAD_non_cancer_OTH_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_ASJ_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_FIN_AF"].value, 0.0),
            math.isclose(maf_record["gnomAD_non_cancer_MAX_AF_adj"].value, 0.000229253),
            maf_record["gnomAD_non_cancer_MAX_AF_POPS_adj"].value == ['amr'],
        ]
    )


def test_variant_with_multiple_max_pop(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_test_vcf_record,
    get_empty_maf_record,
):

    # setup annotator
    ref_path = get_test_file("gnomad")
    ref_pattern = 'gnomad_test.{}.feather'
    annotator = setup_annotator(test_scheme, ref_path, ref_pattern)

    # setup vcf record
    vcf_record = get_test_vcf_record(
        chrom="chr4", pos=10568, stop=10568, ref="T", alleles=("T", "G"), alts=("G",)
    )

    # annotate maf record
    maf_record = annotator.annotate(get_empty_maf_record, vcf_record)

    print(
        [
            maf_record["gnomAD_non_cancer_EAS_AF"].value,
            maf_record["gnomAD_non_cancer_AFR_AF"].value,
            maf_record["gnomAD_non_cancer_AMI_AF"].value,
            maf_record["gnomAD_non_cancer_MID_AF"].value,
            maf_record["gnomAD_non_cancer_SAS_AF"].value,
            maf_record["gnomAD_non_cancer_NFE_AF"].value,
            maf_record["gnomAD_non_cancer_AF"].value,
            maf_record["gnomAD_non_cancer_AMR_AF"].value,
            maf_record["gnomAD_non_cancer_OTH_AF"].value,
            maf_record["gnomAD_non_cancer_ASJ_AF"].value,
            maf_record["gnomAD_non_cancer_FIN_AF"].value,
            maf_record["gnomAD_non_cancer_MAX_AF_adj"].value,
            maf_record["gnomAD_non_cancer_MAX_AF_POPS_adj"].value,
        ]
    )

    assert all(
        [
            math.isclose(maf_record["gnomAD_non_cancer_EAS_AF"].value, 0.999461),
            math.isclose(maf_record["gnomAD_non_cancer_AFR_AF"].value, 0.998853),
            math.isclose(maf_record["gnomAD_non_cancer_AMI_AF"].value, 1),
            math.isclose(maf_record["gnomAD_non_cancer_MID_AF"].value, 1),
            math.isclose(maf_record["gnomAD_non_cancer_SAS_AF"].value, 0.997885),
            math.isclose(maf_record["gnomAD_non_cancer_NFE_AF"].value, 0.999238),
            math.isclose(maf_record["gnomAD_non_cancer_AF"].value, 0.999032),
            math.isclose(maf_record["gnomAD_non_cancer_AMR_AF"].value, 0.999063),
            math.isclose(maf_record["gnomAD_non_cancer_OTH_AF"].value, 0.998087),
            math.isclose(maf_record["gnomAD_non_cancer_ASJ_AF"].value, 0.999326),
            math.isclose(maf_record["gnomAD_non_cancer_FIN_AF"].value, 0.998418),
            math.isclose(maf_record["gnomAD_non_cancer_MAX_AF_adj"].value, 1),
            maf_record["gnomAD_non_cancer_MAX_AF_POPS_adj"].value == ['ami', 'mid'],
        ]
    )


def test_variant_on_nonexistant_chromosome(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_test_vcf_record,
    get_empty_maf_record,
):

    # setup annotator
    ref_path = get_test_file("gnomad")
    ref_pattern = 'gnomad_test.{}.feather'
    annotator = setup_annotator(test_scheme, ref_path, ref_pattern)

    # setup vcf record
    vcf_record = get_test_vcf_record(
        chrom="I don't exist",
        pos=10568,
        stop=10568,
        ref="T",
        alleles=("T", "G"),
        alts=("G",),
    )

    # annotate maf record
    with pytest.raises(KeyError, match=r'Unrecognized contig encountered: .+'):
        maf_record = annotator.annotate(get_empty_maf_record, vcf_record)
