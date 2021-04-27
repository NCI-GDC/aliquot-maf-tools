"""
Tests for the ``aliquotmaf.filters.OffTarget`` class.
"""
from collections import OrderedDict

import pytest
from maflib.column_types import OneBasedIntegerColumn, StringColumn

from aliquotmaf.converters.builder import get_builder
from aliquotmaf.filters import OffTarget


@pytest.fixture
def setup_filter():
    created = []

    def _make_filter(source):
        curr = OffTarget.setup(source)
        created.append(curr)
        return curr

    yield _make_filter

    for record in created:
        record.shutdown()


@pytest.fixture
def test_scheme(get_test_scheme):
    vals = [("vcf_region", StringColumn), ("End_Position", OneBasedIntegerColumn)]

    coldict = OrderedDict(vals)
    return get_test_scheme(coldict)


def test_setup_offtarget(setup_filter, get_test_file):
    bed_file = get_test_file("fake_regions.bed.gz")
    filterer = setup_filter([bed_file])
    assert isinstance(filterer, OffTarget)


@pytest.mark.parametrize(
    "vcf_region, end_position, expected",
    [
        ("chr1:11:.:.:.", 11, True),
        ("chr1:11800:.:.:.", 11868, True),
        ("chr1:11800:.:.:.", 11869, False),
        ("chr1:12227:.:.:.", 12228, False),
        ("chr2:11800:.:.:.", 11869, True),
        ("chr1:939039:.:.:.", 939040, False),
        ("chr10:945517:.:.:.", 945520, False),
    ],
)
def test_offtarget_filter(
    test_scheme,
    setup_filter,
    get_test_file,
    get_empty_maf_record,
    vcf_region,
    end_position,
    expected,
):
    """
    Test offtarget filter 
    """
    bed_file = [
        get_test_file("fake_regions.bed.gz"),
        get_test_file("fake_regions_2.bed.gz"),
    ]
    filterer = setup_filter(bed_file)
    maf_record = get_empty_maf_record
    maf_record["vcf_region"] = get_builder("vcf_region", test_scheme, value=vcf_region)
    maf_record["End_Position"] = get_builder(
        "End_Position", test_scheme, value=end_position
    )
    result = filterer.filter(maf_record)
    assert result is expected
