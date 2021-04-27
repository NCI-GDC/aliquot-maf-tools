"""
Tests for the ``aliquotmaf.annotators.DbSnpValidation`` class.
"""
from collections import OrderedDict

import pysam
import pytest
from maflib.column_types import SequenceOfStrings
from maflib.record import MafColumnRecord

from aliquotmaf.annotators import DbSnpValidation
from aliquotmaf.converters.builder import get_builder


@pytest.fixture
def setup_annotator():
    created = []

    def _make_annotator(scheme, source):
        curr = DbSnpValidation.setup(scheme, source)
        created.append(curr)
        return curr

    yield _make_annotator

    for record in created:
        record.shutdown()


@pytest.fixture
def test_scheme(get_test_scheme):
    coldict = OrderedDict(
        [("dbSNP_RS", SequenceOfStrings), ("dbSNP_Val_Status", SequenceOfStrings)]
    )
    return get_test_scheme(coldict)


def test_setup_dbsnp(test_scheme, setup_annotator, get_test_file):
    db_path = get_test_file("dbsnp_valstatus.db")
    annotator = setup_annotator(test_scheme, source=db_path)
    assert isinstance(annotator, DbSnpValidation)


def test_annotate_dbsnp(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record
):
    db_path = get_test_file("dbsnp_valstatus.db")
    annotator = setup_annotator(test_scheme, source=db_path)

    ## should match
    maf_record = get_empty_maf_record
    maf_record["dbSNP_RS"] = get_builder("dbSNP_RS", test_scheme, value="rs540")
    maf_record = annotator.annotate(maf_record)

    assert maf_record["dbSNP_Val_Status"].value == ["byOtherPop"]

    ## novel should return empty list
    maf_record["dbSNP_RS"] = get_builder("dbSNP_RS", test_scheme, value="novel")
    maf_record["dbSNP_Val_Status"] = get_builder(
        "dbSNP_Val_Status", test_scheme, value=None
    )
    maf_record = annotator.annotate(maf_record)

    assert maf_record["dbSNP_Val_Status"].value == []

    ## empty should return empty list
    maf_record["dbSNP_RS"] = get_builder("dbSNP_RS", test_scheme, value=None)
    maf_record["dbSNP_Val_Status"] = get_builder(
        "dbSNP_Val_Status", test_scheme, value=None
    )
    maf_record = annotator.annotate(maf_record)

    assert maf_record["dbSNP_Val_Status"].value == []
