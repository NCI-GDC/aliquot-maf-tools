"""
Tests for the ``aliquotmaf.filters.FilterGnomAD`` class.
"""

from collections import OrderedDict

import pytest
from maflib.column_types import NullableFloatColumn, SequenceOfStrings

from aliquotmaf.converters.builder import get_builder
from aliquotmaf.filters import FilterGnomAD

GNOMAD_MAF_COLUMNS = (
    "gnomAD_non_cancer_EAS_AF",
    "gnomAD_non_cancer_AFR_AF",
    "gnomAD_non_cancer_AMI_AF",
    "gnomAD_non_cancer_MID_AF",
    "gnomAD_non_cancer_SAS_AF",
    "gnomAD_non_cancer_NFE_AF",
    "gnomAD_non_cancer_AF",
    "gnomAD_non_cancer_AMR_AF",
    "gnomAD_non_cancer_OTH_AF",
    "gnomAD_non_cancer_ASJ_AF",
    "gnomAD_non_cancer_FIN_AF",
    "gnomAD_non_cancer_MAX_AF_adj",
    "gnomAD_non_cancer_MAX_AF_POPS_adj",
)


@pytest.fixture
def setup_filter():
    created = []

    def _make_filter(cutoff):
        curr = FilterGnomAD.setup(cutoff)
        created.append(curr)
        return curr

    yield _make_filter

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


def test_setup_filter_gnomAD(setup_filter):
    cutoff = 0.0004
    filterer = setup_filter(cutoff)
    assert isinstance(filterer, FilterGnomAD)


def test_filter_gnomad_null(test_scheme, setup_filter, get_empty_maf_record):
    """
    Test FilterGnomAD when data is null
    """
    cutoff = 0.0004
    filterer = setup_filter(cutoff)
    maf_record = get_empty_maf_record
    for key in GNOMAD_MAF_COLUMNS:
        value = ""
        if key == "gnomAD_non_cancer_MAX_AF_POPS_adj":
            value = []
        maf_record[key] = get_builder(key, test_scheme, value=value)
    result = filterer.filter(maf_record)
    assert result is False


def test_filter_gnomad_positive(test_scheme, setup_filter, get_empty_maf_record):
    """
    Test FilterGnomAD when value is above cutoff
    """
    cutoff = 0.0004
    filterer = setup_filter(cutoff)
    maf_record = get_empty_maf_record
    for key in GNOMAD_MAF_COLUMNS:
        value = ""
        if key == "gnomAD_non_cancer_MAX_AF_adj":
            value = 0.5
        if key == "gnomAD_non_cancer_MAX_AF_POPS_adj":
            value = []
        maf_record[key] = get_builder(key, test_scheme, value=value)
    result = filterer.filter(maf_record)
    assert result is True


def test_filter_gnomad_on_cutoff(test_scheme, setup_filter, get_empty_maf_record):
    """
    Test FilterGnomAD when value is exactly cutoff
    """
    cutoff = 0.0004
    filterer = setup_filter(cutoff)
    maf_record = get_empty_maf_record
    for key in GNOMAD_MAF_COLUMNS:
        value = ""
        if key == "gnomAD_non_cancer_MAX_AF_adj":
            value = cutoff
        if key == "gnomAD_non_cancer_MAX_AF_POPS_adj":
            value = []
        maf_record[key] = get_builder(key, test_scheme, value=value)
    result = filterer.filter(maf_record)
    assert result is False


def test_filter_gnomad_negative(test_scheme, setup_filter, get_empty_maf_record):
    """
    Test FilterGnomAD when data is below cutoff
    """
    cutoff = 0.0004
    filterer = setup_filter(cutoff)
    maf_record = get_empty_maf_record
    for key in GNOMAD_MAF_COLUMNS:
        value = ""
        if key == "gnomAD_non_cancer_MAX_AF_adj":
            value = 0.000399
        if key == "gnomAD_non_cancer_MAX_AF_POPS_adj":
            value = []
        maf_record[key] = get_builder(key, test_scheme, value=value)
    result = filterer.filter(maf_record)
    assert result is False
