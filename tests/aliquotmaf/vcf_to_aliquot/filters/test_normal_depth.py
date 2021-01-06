#!/usr/bin/env python3
"""
Tests for the ``aliquotmaf.filters.NormalDepth`` class.
"""
from collections import OrderedDict
from types import SimpleNamespace

import pytest
from maflib.column_types import NullableZeroBasedIntegerColumn

from aliquotmaf.vcf_to_aliquot.converters.builder import get_builder
from aliquotmaf.vcf_to_aliquot.filters import normal_depth as MOD


@pytest.fixture
def setup_filter():
    created = []

    def _make_filter(args):
        curr = MOD.NormalDepth.setup(args)
        created.append(curr)
        return curr

    yield _make_filter

    for record in created:
        record.shutdown()


@pytest.fixture
def test_scheme(get_test_scheme):
    vals = [("n_depth", NullableZeroBasedIntegerColumn)]

    coldict = OrderedDict(vals)
    return get_test_scheme(coldict)


def test_setup_normal_depth(setup_filter):
    args = SimpleNamespace(min_n_depth=7, is_tumor_only=False)
    filterer = setup_filter(args)
    assert isinstance(filterer, MOD.NormalDepth)


@pytest.mark.parametrize(
    "normal_depth, expected", [(None, False), (7, True), (8, False),],
)
def test_normal_depth_filter(
    test_scheme, setup_filter, get_empty_maf_record, normal_depth, expected
):
    """
    Test Normal Depth filter
    """
    args = SimpleNamespace(min_n_depth=7, is_tumor_only=False)
    filterer = setup_filter(args)
    maf_record = get_empty_maf_record
    maf_record["n_depth"] = get_builder("n_depth", test_scheme, value=normal_depth)
    result = filterer.filter(maf_record)
    assert result is expected


# __END__
