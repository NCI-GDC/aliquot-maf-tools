"""
Tests for the ``aliquotmaf.annotators.Hotpot`` class.
"""

from collections import OrderedDict

import pytest
from maflib.column_types import NullableStringColumn, NullableYOrN, StringColumn

from aliquotmaf.annotators import Hotspot
from aliquotmaf.converters.builder import get_builder


@pytest.fixture
def setup_annotator():
    created = []

    def _make_annotator(scheme, source):
        curr = Hotspot.setup(scheme, source)
        created.append(curr)
        return curr

    yield _make_annotator

    for record in created:
        record.shutdown()


@pytest.fixture
def test_scheme(get_test_scheme):
    vals = [
        ("Hugo_Symbol", StringColumn),
        ("HGVSp_Short", NullableStringColumn),
        ("hotspot", NullableYOrN),
    ]

    coldict = OrderedDict(vals)
    return get_test_scheme(coldict)


def test_setup_hotspot(test_scheme, setup_annotator, get_test_file):
    tsv_path = get_test_file("fake_hotspot.tsv")
    annotator = setup_annotator(test_scheme, source=tsv_path)
    assert isinstance(annotator, Hotspot)


def test_hotspot_annotator_no_overlap_1(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record
):
    """
    Not a hotspot
    """
    tsv_path = get_test_file("fake_hotspot.tsv")
    annotator = setup_annotator(test_scheme, source=tsv_path)

    maf_record = get_empty_maf_record
    maf_record["Hugo_Symbol"] = get_builder("Hugo_Symbol", test_scheme, value="Unknown")
    maf_record["HGVSp_Short"] = get_builder("HGVSp_Short", test_scheme, value=None)

    maf_record = annotator.annotate(maf_record)

    assert maf_record["hotspot"].value.value == "N"


def test_hotspot_annotator_no_overlap_2(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record
):
    """
    Not a hotspot
    """
    tsv_path = get_test_file("fake_hotspot.tsv")
    annotator = setup_annotator(test_scheme, source=tsv_path)

    maf_record = get_empty_maf_record
    maf_record["Hugo_Symbol"] = get_builder("Hugo_Symbol", test_scheme, value="ASXL1")
    maf_record["HGVSp_Short"] = get_builder(
        "HGVSp_Short", test_scheme, value="p.R547fs"
    )

    maf_record = annotator.annotate(maf_record)

    assert maf_record["hotspot"].value.value == "N"


def test_hotspot_annotator_overlap(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record
):
    """
    Is a hotspot
    """
    tsv_path = get_test_file("fake_hotspot.tsv")
    annotator = setup_annotator(test_scheme, source=tsv_path)

    maf_record = get_empty_maf_record
    maf_record["Hugo_Symbol"] = get_builder("Hugo_Symbol", test_scheme, value="ASXL1")
    maf_record["HGVSp_Short"] = get_builder(
        "HGVSp_Short", test_scheme, value="p.R548fs"
    )

    maf_record = annotator.annotate(maf_record)

    assert maf_record["hotspot"].value.value == "Y"
