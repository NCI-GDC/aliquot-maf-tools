"""
Tests for the ``aliquotmaf.filters.ExAC`` class.
"""
import pytest
from collections import OrderedDict

from maflib.column_types import NullableFloatColumn

from aliquotmaf.filters import ExAC
from aliquotmaf.converters.builder import get_builder

subpops = [
    "nontcga_ExAC_AF_Adj",
    "nontcga_ExAC_AF",
    "nontcga_ExAC_AF_AFR",
    "nontcga_ExAC_AF_AMR",
    "nontcga_ExAC_AF_EAS",
    "nontcga_ExAC_AF_FIN",
    "nontcga_ExAC_AF_NFE",
    "nontcga_ExAC_AF_OTH",
    "nontcga_ExAC_AF_SAS",
]


@pytest.fixture
def setup_filter():
    created = []

    def _make_filter(cutoff):
        curr = ExAC.setup(cutoff)
        created.append(curr)
        return curr

    yield _make_filter

    for record in created:
        record.shutdown()


@pytest.fixture
def test_scheme(get_test_scheme):
    vals = []
    for p in subpops:
        vals.append((p, NullableFloatColumn))

    coldict = OrderedDict(vals)
    return get_test_scheme(coldict)


def test_setup_exac(setup_filter):
    cutoff = 0.0004
    filterer = setup_filter(cutoff)
    assert isinstance(filterer, ExAC)


def test_exac_filter_1(test_scheme, setup_filter, get_empty_maf_record):
    """
    Test exac filter when all freqs are None
    """
    cutoff = 0.0004
    filterer = setup_filter(cutoff)
    maf_record = get_empty_maf_record
    for key in subpops:
        maf_record[key] = get_builder(key, test_scheme, value=None)
    result = filterer.filter(maf_record)
    assert result is False


def test_exac_filter_2(test_scheme, setup_filter, get_empty_maf_record):
    """
    Test exac filter when all freqs are below cutoff 
    """
    cutoff = 0.0004
    filterer = setup_filter(cutoff)
    maf_record = get_empty_maf_record
    for key in subpops:
        maf_record[key] = get_builder(key, test_scheme, value=0.0003)
    result = filterer.filter(maf_record)
    assert result is False


def test_exac_filter_3(test_scheme, setup_filter, get_empty_maf_record):
    """
    Test exac filter when all freqs are exactly cutoff 
    """
    cutoff = 0.0004
    filterer = setup_filter(cutoff)
    maf_record = get_empty_maf_record
    for key in subpops:
        maf_record[key] = get_builder(key, test_scheme, value=0.0004)
    result = filterer.filter(maf_record)
    assert result is False


def test_exac_filter_4(test_scheme, setup_filter, get_empty_maf_record):
    """
    Test exac filter when all but 1 freqs is above cutoff 
    """
    cutoff = 0.0004
    filterer = setup_filter(cutoff)
    maf_record = get_empty_maf_record
    for key in subpops:
        maf_record[key] = get_builder(key, test_scheme, value=0.0004)
    maf_record["nontcga_ExAC_AF_Adj"] = get_builder(
        "nontcga_ExAC_AF_Adj", test_scheme, value=0.00041
    )
    result = filterer.filter(maf_record)
    assert result is True
