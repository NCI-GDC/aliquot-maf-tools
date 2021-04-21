"""
Tests for the ``aliquotmaf.filters.GdcPon`` class.
"""
from collections import OrderedDict

import pytest
from maflib.column_types import StringColumn

from aliquotmaf.converters.builder import get_builder
from aliquotmaf.filters import GdcPon


@pytest.fixture
def setup_filter():
    created = []

    def _make_filter(source):
        curr = GdcPon.setup(source)
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


def test_setup_pon(setup_filter, get_test_file):
    vcf_path = get_test_file("fake_exac.vcf.gz")
    filterer = setup_filter(vcf_path)
    assert isinstance(filterer, GdcPon)


@pytest.mark.parametrize(
    "vcf_region, expected",
    [
        ("chr1:11:.:G:C", False),
        ("chr1:10:.:C:T", True),
        ("chr2:8:.:CTACTT:C", False),
        ("chr2:10:.:A:C", True),
    ],
)
def test_pon_filter(
    test_scheme, setup_filter, get_test_file, get_empty_maf_record, vcf_region, expected
):
    """
    Test pon filter
    """
    vcf_path = get_test_file("fake_exac.vcf.gz")
    filterer = setup_filter(vcf_path)
    maf_record = get_empty_maf_record
    maf_record["vcf_region"] = get_builder("vcf_region", test_scheme, value=vcf_region)
    result = filterer.filter(maf_record)
    assert result is expected
