"""
Tests for the ``aliquotmaf.annotators.Cosmic`` class.
"""
from collections import OrderedDict, namedtuple

import pytest
from maflib.column_types import EntrezGeneId, NullableStringColumn, StringColumn
from maflib.record import MafColumnRecord

from aliquotmaf.annotators import MAF_FEATURE, MAF_SYMBOL, Entrez
from aliquotmaf.converters.builder import get_builder


@pytest.fixture
def setup_annotator():
    created = []

    def _make_annotator(scheme, entrez_json_file):
        curr = Entrez.setup(scheme, entrez_json_file)
        created.append(curr)
        return curr

    yield _make_annotator

    for record in created:
        record.shutdown()


@pytest.fixture
def test_scheme(get_test_scheme):
    coldict = OrderedDict(
        [
            (MAF_SYMBOL, StringColumn),
            (MAF_FEATURE, NullableStringColumn),
            ("Entrez_Gene_Id", EntrezGeneId),
        ]
    )
    return get_test_scheme(coldict)


def test_setup_entrez(test_scheme, setup_annotator, get_test_file):
    json_path = get_test_file("ex_entrez.json")
    annotator = setup_annotator(test_scheme, entrez_json_file=json_path)
    assert isinstance(annotator, Entrez)


def test_entrez_symbol_and_feature(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record,
):

    # setup annotator
    json_path = get_test_file("ex_entrez.json")
    annotator = setup_annotator(test_scheme, entrez_json_file=json_path)

    init_maf_record = get_empty_maf_record
    init_maf_record[MAF_SYMBOL] = get_builder(
        MAF_SYMBOL, test_scheme, value='PRAMEF27', default=''
    )
    init_maf_record[MAF_FEATURE] = get_builder(
        MAF_FEATURE, test_scheme, value='ENST00000436041', default=''
    )

    # print(test_scheme.column_class('Entrez_Gene_Id').__name__)
    maf_record = annotator.annotate(init_maf_record)

    assert maf_record['Entrez_Gene_Id'].value == 101929983


def test_entrez_symbol_only(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record,
):

    # setup annotator
    json_path = get_test_file("ex_entrez.json")
    annotator = setup_annotator(test_scheme, entrez_json_file=json_path)

    init_maf_record = get_empty_maf_record
    init_maf_record[MAF_SYMBOL] = get_builder(
        MAF_SYMBOL, test_scheme, value='PRAMEF27', default=''
    )
    # print(test_scheme.column_class('Entrez_Gene_Id').__name__)
    maf_record = annotator.annotate(init_maf_record)

    assert maf_record['Entrez_Gene_Id'].value == 101929983


def test_feature_only(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record,
):

    # setup annotator
    json_path = get_test_file("ex_entrez.json")
    annotator = setup_annotator(test_scheme, entrez_json_file=json_path)

    init_maf_record = get_empty_maf_record

    init_maf_record[MAF_FEATURE] = get_builder(
        MAF_FEATURE, test_scheme, value='ENST00000436041', default=''
    )

    # print(test_scheme.column_class('Entrez_Gene_Id').__name__)
    maf_record = annotator.annotate(init_maf_record)

    assert maf_record['Entrez_Gene_Id'].value == 101929983


def test_neither_id_present_in_query(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record,
):

    # setup annotator
    json_path = get_test_file("ex_entrez.json")
    annotator = setup_annotator(test_scheme, entrez_json_file=json_path)

    init_maf_record = get_empty_maf_record

    # print(test_scheme.column_class('Entrez_Gene_Id').__name__)
    maf_record = annotator.annotate(init_maf_record)

    assert maf_record['Entrez_Gene_Id'].value is None


def test_symbol_not_present_in_database(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record,
):
    # setup annotator
    json_path = get_test_file("ex_entrez.json")
    annotator = setup_annotator(test_scheme, entrez_json_file=json_path)

    init_maf_record = get_empty_maf_record
    init_maf_record[MAF_SYMBOL] = get_builder(
        MAF_SYMBOL, test_scheme, value='NOTAGENE', default=''
    )

    # print(test_scheme.column_class('Entrez_Gene_Id').__name__)
    maf_record = annotator.annotate(init_maf_record)

    assert maf_record['Entrez_Gene_Id'].value is None


def test_gencode_id_not_present_in_database(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record,
):
    # setup annotator
    json_path = get_test_file("ex_entrez.json")
    annotator = setup_annotator(test_scheme, entrez_json_file=json_path)

    init_maf_record = get_empty_maf_record
    init_maf_record[MAF_FEATURE] = get_builder(
        MAF_FEATURE, test_scheme, value='ENST99999999999', default=''
    )

    # print(test_scheme.column_class('Entrez_Gene_Id').__name__)
    maf_record = annotator.annotate(init_maf_record)

    assert maf_record['Entrez_Gene_Id'].value is None
