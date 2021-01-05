#!/usr/bin/env python3
"""
Tests for the ``aliquotmaf.annotators.Hotpot`` class.
"""
from collections import OrderedDict
from types import SimpleNamespace

import pytest
from maflib.column_types import NullableStringColumn, NullableYOrN, StringColumn
from maflib.record import MafColumnRecord

from aliquotmaf.vcf_to_aliquot.annotators import hotspot as MOD
from aliquotmaf.vcf_to_aliquot.converters.builder import get_builder


@pytest.fixture
def setup_annotator():
    created = []

    def _make_annotator(scheme, args):
        curr = MOD.Hotspot.setup(scheme, args)
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
    args = SimpleNamespace(hotspot_tsv=tsv_path)
    annotator = setup_annotator(test_scheme, args)
    assert isinstance(annotator, MOD.Hotspot)


def test_hotspot_annotator_no_overlap_1(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record
):
    """
    Not a hotspot 
    """
    tsv_path = get_test_file("fake_hotspot.tsv")
    args = SimpleNamespace(hotspot_tsv=tsv_path)
    annotator = setup_annotator(test_scheme, args)

    maf_record = get_empty_maf_record
    maf_record["Hugo_Symbol"] = get_builder("Hugo_Symbol", test_scheme, value="Unknown")
    maf_record["HGVSp_Short"] = get_builder("HGVSp_Short", test_scheme, value=None)

    found = annotator.annotate(maf_record)

    assert found.value.value == "N"


def test_hotspot_annotator_no_overlap_2(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record
):
    """
    Not a hotspot 
    """
    tsv_path = get_test_file("fake_hotspot.tsv")
    args = SimpleNamespace(hotspot_tsv=tsv_path)
    annotator = setup_annotator(test_scheme, args)

    maf_record = get_empty_maf_record
    maf_record["Hugo_Symbol"] = get_builder("Hugo_Symbol", test_scheme, value="ASXL1")
    maf_record["HGVSp_Short"] = get_builder(
        "HGVSp_Short", test_scheme, value="p.R547fs"
    )

    found = annotator.annotate(maf_record)

    assert found.value.value == "N"


def test_hotspot_annotator_overlap(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record
):
    """
    Is a hotspot 
    """
    tsv_path = get_test_file("fake_hotspot.tsv")
    args = SimpleNamespace(hotspot_tsv=tsv_path)
    annotator = setup_annotator(test_scheme, args)

    maf_record = get_empty_maf_record
    maf_record["Hugo_Symbol"] = get_builder("Hugo_Symbol", test_scheme, value="ASXL1")
    maf_record["HGVSp_Short"] = get_builder(
        "HGVSp_Short", test_scheme, value="p.R548fs"
    )

    found = annotator.annotate(maf_record)

    assert found.value.value == "Y"


# __END__
