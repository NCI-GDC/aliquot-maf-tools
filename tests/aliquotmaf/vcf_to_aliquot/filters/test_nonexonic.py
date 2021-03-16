#!/usr/bin/env python3
"""
Tests for the ``aliquotmaf.filters.NonExonic`` class.
"""
from collections import OrderedDict
from types import SimpleNamespace

import pytest
from maflib.column_types import OneBasedIntegerColumn, StringColumn

from aliquotmaf.vcf_to_aliquot.converters.builder import get_builder
from aliquotmaf.vcf_to_aliquot.filters import nonexonic as MOD


@pytest.fixture
def setup_filter():
    created = []

    def _make_filter(args):
        curr = MOD.NonExonic.setup(args)
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


def test_setup_nonexonic(setup_filter, get_test_file):
    bed_file = get_test_file("fake_regions.bed.gz")
    args = SimpleNamespace(nonexonic_intervals=bed_file)
    filterer = setup_filter(args)
    assert isinstance(filterer, MOD.NonExonic)


@pytest.mark.parametrize(
    "vcf_region, end_position, expected",
    [
        ("chr1:11:.:.:.", 11, True),
        ("chr1:11800:.:.:.", 11868, True),
        ("chr1:11800:.:.:.", 11869, False),
        ("chr1:12227:.:.:.", 12228, False),
        ("chr2:11800:.:.:.", 11869, True),
    ],
)
def test_nonexonic_filter(
    test_scheme,
    setup_filter,
    get_test_file,
    get_empty_maf_record,
    vcf_region,
    end_position,
    expected,
):
    """
    Test nonexonic filter
    """
    bed_file = get_test_file("fake_regions.bed.gz")
    args = SimpleNamespace(nonexonic_intervals=bed_file)
    filterer = setup_filter(args)
    maf_record = get_empty_maf_record
    maf_record["vcf_region"] = get_builder("vcf_region", test_scheme, value=vcf_region)
    maf_record["End_Position"] = get_builder(
        "End_Position", test_scheme, value=end_position
    )
    result = filterer.filter(maf_record)
    assert result is expected


# __END__
