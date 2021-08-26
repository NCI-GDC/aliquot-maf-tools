"""
Tests for the ``aliquotmaf.filters.NormalDepth`` class.
"""
from collections import OrderedDict

import pytest
from maflib.column_types import NullableZeroBasedIntegerColumn

from aliquotmaf.converters.builder import get_builder
from aliquotmaf.filters import NormalDepth


@pytest.fixture
def setup_filter():
    created = []

    def _make_filter(cutoff):
        curr = NormalDepth.setup(cutoff)
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
    filterer = setup_filter(7)
    assert isinstance(filterer, NormalDepth)


@pytest.mark.parametrize(
    "normal_depth, expected", [(None, False), (7, True), (8, False)]
)
def test_normal_depth_filter(
    test_scheme, setup_filter, get_empty_maf_record, normal_depth, expected
):
    """
    Test Normal Depth filter
    """
    filterer = setup_filter(7)
    maf_record = get_empty_maf_record
    maf_record["n_depth"] = get_builder("n_depth", test_scheme, value=normal_depth)
    result = filterer.filter(maf_record)
    assert result is expected
