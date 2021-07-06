"""
Tests for the ``aliquotmaf.annotators.Cosmic`` class.
"""
from collections import OrderedDict, namedtuple

import pytest
from maflib.column_types import EntrezGeneId
from maflib.record import MafColumnRecord

from aliquotmaf.annotators import Entrez, Vcf
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
    coldict = OrderedDict([("Entrez_Gene_Id", EntrezGeneId)])
    return get_test_scheme(coldict)


def test_setup_entrez(test_scheme, setup_annotator, get_test_file):
    json_path = get_test_file("ex_entrez.json")
    annotator = setup_annotator(test_scheme, entrez_json_file=json_path)
    assert isinstance(annotator, Entrez)


def test_entrez_symbol_and_feature(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_test_vcf_record,
    get_empty_maf_record,
):

    # setup annotator
    json_path = get_test_file("ex_entrez.json")
    annotator = setup_annotator(test_scheme, entrez_json_file=json_path)

    vcf_record = get_test_vcf_record(
        chrom="chr1",
        pos=13053533,
        stop=13053533,
        ref="T",
        alleles=("T", "G"),
        alts=("G",),
        info={Vcf.SYMBOL: 'PRAMEF27', Vcf.Feature: 'ENST00000436041'},
    )
    # print(test_scheme.column_class('Entrez_Gene_Id').__name__)
    maf_record = annotator.annotate(get_empty_maf_record, vcf_record)

    assert maf_record['Entrez_Gene_Id'].value == 101929983


def test_entrez_symbol(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_test_vcf_record,
    get_empty_maf_record,
):

    # setup annotator
    json_path = get_test_file("ex_entrez.json")
    annotator = setup_annotator(test_scheme, entrez_json_file=json_path)

    vcf_record = get_test_vcf_record(
        chrom="chr1",
        pos=13053533,
        stop=13053533,
        ref="T",
        alleles=("T", "G"),
        alts=("G",),
        info={Vcf.SYMBOL: 'PRAMEF27', Vcf.Feature: 'ENST00000436041'},
    )
    # print(test_scheme.column_class('Entrez_Gene_Id').__name__)
    maf_record = annotator.annotate(get_empty_maf_record, vcf_record)

    assert maf_record['Entrez_Gene_Id'].value == 101929983
