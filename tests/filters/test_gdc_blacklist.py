"""
Tests for the ``aliquotmaf.filters.GdcBlacklist`` class.
"""

from collections import OrderedDict

import pytest
from maflib.column_types import UUIDColumn

from aliquotmaf.converters.builder import get_builder
from aliquotmaf.filters import GdcBlacklist


@pytest.fixture
def setup_filter():
    created = []

    def _make_filter(source):
        curr = GdcBlacklist.setup(source)
        created.append(curr)
        return curr

    yield _make_filter

    for record in created:
        record.shutdown()


@pytest.fixture
def test_scheme(get_test_scheme):
    vals = [("Tumor_Sample_UUID", UUIDColumn)]

    coldict = OrderedDict(vals)
    return get_test_scheme(coldict)


def test_setup_blacklist(setup_filter, get_test_file):
    tsv_path = get_test_file("fake_blacklist.tsv")
    filterer = setup_filter(tsv_path)
    assert isinstance(filterer, GdcBlacklist)


@pytest.mark.parametrize(
    "tumor_uuid, expected_bool, expected_tags",
    [
        ("00000000-0000-0000-0000-000000000002", False, []),
        ("00000000-0000-0000-0000-000000000000", True, ["QC_Pending"]),
        ("00000000-0000-0000-0000-000000000001", True, ["QC_Pending", "OTHER"]),
    ],
)
def test_blacklist_filter(
    test_scheme,
    setup_filter,
    get_test_file,
    get_empty_maf_record,
    tumor_uuid,
    expected_bool,
    expected_tags,
):
    """
    Test blacklist filter
    """
    tsv_path = get_test_file("fake_blacklist.tsv")
    filterer = setup_filter(tsv_path)
    maf_record = get_empty_maf_record
    maf_record["Tumor_Sample_UUID"] = get_builder(
        "Tumor_Sample_UUID", test_scheme, value=tumor_uuid
    )
    result = filterer.filter(maf_record)
    assert result is expected_bool
    assert filterer.tags == expected_tags
