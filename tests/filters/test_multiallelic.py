"""
Tests for the ``aliquotmaf.filters.Multiallelic`` class.
"""
import pytest
from collections import OrderedDict

from maflib.column_types import StringColumn

from aliquotmaf.filters import Multiallelic
from aliquotmaf.converters.builder import get_builder


@pytest.fixture
def setup_filter():
    created = []

    def _make_filter():
        curr = Multiallelic.setup()
        created.append(curr)
        return curr

    yield _make_filter

    for record in created:
        record.shutdown()


@pytest.fixture
def test_scheme(get_test_scheme):
    vals = [("vcf_region", StringColumn)]

    coldict = OrderedDict(vals)
    return get_test_scheme(coldict)


def test_setup_multiallelic(setup_filter):
    filterer = setup_filter()
    assert isinstance(filterer, Multiallelic)


@pytest.mark.parametrize(
    "vcf_region, expected",
    [
        ("chr1:11:.:G:C", False),
        ("chr1:10:.:C:T,G", True),
        ("chr1:10:.:C:T,G,A", True),
        ("chr2:8:.:CTACTT:C", False),
        ("chr1:10:.:C:T,T,C", False),
    ],
)
def test_multiallelic_filter(
    test_scheme, setup_filter, get_empty_maf_record, vcf_region, expected
):
    """
    Test multiallelic filter 
    """
    filterer = setup_filter()
    maf_record = get_empty_maf_record
    maf_record["vcf_region"] = get_builder("vcf_region", test_scheme, value=vcf_region)
    result = filterer.filter(maf_record)
    assert result is expected
