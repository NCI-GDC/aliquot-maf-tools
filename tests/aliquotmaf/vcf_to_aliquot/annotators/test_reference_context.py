#!/usr/bin/env python3
"""
Tests for the ``aliquotmaf.annotators.ReferenceContext`` class.
"""
from collections import OrderedDict, namedtuple
from types import SimpleNamespace

import pysam
import pytest
from maflib.column_types import StringColumn

from aliquotmaf.vcf_to_aliquot.annotators import reference_context as MOD


@pytest.fixture
def setup_annotator():
    created = []

    def _make_annotator(scheme, args, context_size=5):
        curr = MOD.ReferenceContext.setup(scheme, args, context_size=context_size)
        created.append(curr)
        return curr

    yield _make_annotator

    for record in created:
        record.shutdown()


@pytest.fixture
def test_scheme(get_test_scheme):
    coldict = OrderedDict([("CONTEXT", StringColumn)])
    return get_test_scheme(coldict)


@pytest.fixture
def vcf_gen(get_test_file):
    VcfRec = namedtuple("VcfRec", ["snp", "deletion", "insertion"])
    created = []

    def _vcf_gen(name):
        vcf_path = get_test_file(name)
        vcf_obj = pysam.VariantFile(vcf_path)
        created.append(vcf_obj)
        gen = vcf_obj.fetch()
        curr = VcfRec(snp=next(gen), deletion=next(gen), insertion=next(gen))
        return curr

    yield _vcf_gen

    for obj in created:
        obj.close()


def test_setup_reference_context(test_scheme, setup_annotator, get_test_file):
    fasta_path = get_test_file("fake_ref.fa")
    args = SimpleNamespace(reference_context=fasta_path)
    annotator = setup_annotator(test_scheme, args)
    assert annotator.context_size == 5


def test_reference_context_snp(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record, vcf_gen
):
    fasta_path = get_test_file("fake_ref.fa")
    args = SimpleNamespace(reference_context=fasta_path)
    annotator = setup_annotator(test_scheme, args)

    gen = vcf_gen("ex1.vcf.gz")

    # simple snp
    record = gen.snp

    ## default context
    maf_record = annotator.annotate(get_empty_maf_record, record)
    assert maf_record["CONTEXT"].value == "AGTGGCTCATT"

    ## truncated left context
    annotator.context_size = 12
    maf_record = annotator.annotate(get_empty_maf_record, record)
    assert maf_record["CONTEXT"].value == "CACTAGTGGCTCATTGTAAATG"

    ## truncated both context
    annotator.context_size = 1575
    maf_record = annotator.annotate(get_empty_maf_record, record)
    assert maf_record["CONTEXT"].value == annotator.fa.fetch(region=record.chrom)


def test_reference_context_del(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record, vcf_gen
):
    fasta_path = get_test_file("fake_ref.fa")
    args = SimpleNamespace(reference_context=fasta_path)
    annotator = setup_annotator(test_scheme, args)

    gen = vcf_gen("ex1.vcf.gz")

    # deletion
    record = gen.deletion
    maf_record = annotator.annotate(get_empty_maf_record, record)
    assert maf_record["CONTEXT"].value == "AATGAACTTCTGTA"


def test_reference_context_ins(
    test_scheme, setup_annotator, get_test_file, get_empty_maf_record, vcf_gen
):
    fasta_path = get_test_file("fake_ref.fa")
    annotator = setup_annotator(test_scheme, source=fasta_path)

    gen = vcf_gen("ex1.vcf.gz")

    # insertion
    record = gen.insertion
    maf_record = annotator.annotate(get_empty_maf_record, record)
    assert maf_record["CONTEXT"].value == "TGTAATTGAAA"


# __END__
