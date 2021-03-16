"""
Tests for the ``aliquotmaf.annotators.DbSnpValidation`` class.
"""
from collections import OrderedDict
from types import SimpleNamespace

import pysam
import pytest
from maflib.column_types import SequenceOfStrings
from maflib.record import MafColumnRecord

from aliquotmaf.vcf_to_aliquot.annotators import dbsnp_validation as MOD
from aliquotmaf.vcf_to_aliquot.converters.builder import get_builder


@pytest.fixture
def setup_annotator():
    created = []

    def _make_annotator(scheme, args):
        curr = MOD.DbSnpValidation.setup(scheme, args)
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
    args = SimpleNamespace(dbsnp_priority_db=db_path)
    annotator = setup_annotator(test_scheme, args)
    assert isinstance(annotator, MOD.DbSnpValidation)


def test_annotate_dbsnp(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record
):
    db_path = get_test_file("dbsnp_valstatus.db")
    args = SimpleNamespace(dbsnp_priority_db=db_path)
    annotator = setup_annotator(test_scheme, args)

    ## should match
    maf_record = get_empty_maf_record
    maf_record["dbSNP_RS"] = get_builder("dbSNP_RS", test_scheme, value="rs540")
    found = annotator.annotate(maf_record)

    assert found.value == ["byOtherPop"]

    ## novel should return empty list
    maf_record["dbSNP_RS"] = get_builder("dbSNP_RS", test_scheme, value="novel")
    maf_record["dbSNP_Val_Status"] = get_builder(
        "dbSNP_Val_Status", test_scheme, value=None
    )
    found = annotator.annotate(maf_record)

    assert found.value == []

    ## empty should return empty list
    maf_record["dbSNP_RS"] = get_builder("dbSNP_RS", test_scheme, value=None)
    maf_record["dbSNP_Val_Status"] = get_builder(
        "dbSNP_Val_Status", test_scheme, value=None
    )
    found = annotator.annotate(maf_record)

    assert found.value == []


# __END__
