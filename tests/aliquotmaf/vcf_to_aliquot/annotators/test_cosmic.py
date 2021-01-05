"""
Tests for the ``aliquotmaf.annotators.Cosmic`` class.
"""
from collections import OrderedDict, namedtuple
from types import SimpleNamespace

import pysam
import pytest
from maflib.column_types import SequenceOfStrings
from maflib.record import MafColumnRecord

from aliquotmaf.vcf_to_aliquot.annotators import cosmic as MOD
from aliquotmaf.vcf_to_aliquot.converters.builder import get_builder


@pytest.fixture
def setup_annotator():
    created = []

    def _make_annotator(scheme, args):
        curr = MOD.CosmicId.setup(scheme, args)
        created.append(curr)
        return curr

    yield _make_annotator

    for record in created:
        record.shutdown()


@pytest.fixture
def test_scheme(get_test_scheme):
    coldict = OrderedDict(
        [("COSMIC", SequenceOfStrings), ("dbSNP_RS", SequenceOfStrings)]
    )
    return get_test_scheme(coldict)


@pytest.fixture
def vcf_gen(get_test_file):
    VcfRec = namedtuple("VcfRec", ["snp1", "deletion", "insertion", "snp2"])
    created = []

    def _vcf_gen(name):
        vcf_path = get_test_file(name)
        vcf_obj = pysam.VariantFile(vcf_path)
        created.append(vcf_obj)
        gen = vcf_obj.fetch()
        curr = VcfRec(
            snp1=next(gen), deletion=next(gen), insertion=next(gen), snp2=next(gen)
        )
        return curr

    yield _vcf_gen

    for obj in created:
        obj.close()


def test_setup_cosmic(test_scheme, setup_annotator, get_test_file):
    vcf_path = get_test_file("ex2.vcf.gz")
    # Mimic CLI arg namespace
    args = SimpleNamespace(cosmic_vcf=vcf_path)
    annotator = setup_annotator(test_scheme, args)
    assert isinstance(annotator, MOD.CosmicId)


def test_cosmic_snp(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_empty_maf_record,
    vcf_gen,
    get_test_vcf_record,
):
    vcf_path = get_test_file("ex2.vcf.gz")
    args = SimpleNamespace(cosmic_vcf=vcf_path)
    annotator = setup_annotator(test_scheme, args)

    gen = vcf_gen("ex2.vcf.gz")

    # simple snp
    record = gen.snp1

    ## setup
    vcf_record = get_test_vcf_record(
        chrom=record.chrom,
        pos=record.pos,
        alleles=record.alleles,
        ref=record.ref,
        alts=record.alts,
    )
    maf_record = get_empty_maf_record
    # maf_record["dbSNP_RS"] = get_builder("dbSNP_RS", test_scheme, value="novel")
    found = annotator.annotate(maf_record, vcf_record, var_allele_idx=1)

    assert found.value == ["COSM0000"]
    # assert return_nt.dbsnp.value == []

    maf_record = get_empty_maf_record
    # maf_record["dbSNP_RS"] = get_builder("dbSNP_RS", test_scheme, value=None)
    found = annotator.annotate(maf_record, vcf_record, var_allele_idx=1)

    assert found.value == ["COSM0000"]
    # assert return_nt.dbsnp.value == []

    record = gen.snp2

    ## setup
    vcf_record = get_test_vcf_record(
        chrom=record.chrom,
        pos=record.pos,
        alleles=record.alleles,
        ref=record.ref,
        alts=record.alts,
    )
    maf_record = get_empty_maf_record
    # maf_record["dbSNP_RS"] = get_builder("dbSNP_RS", test_scheme, value="novel")
    found = annotator.annotate(maf_record, vcf_record, var_allele_idx=1)

    # Do overwrite dnSNP_RS with cosmic IDs
    assert found.value == ["COSM0003"]
    # assert return_nt.dbsnp.value == []

    ## setup no overlap alleles
    vcf_record = get_test_vcf_record(
        chrom=record.chrom,
        pos=record.pos,
        alleles=(record.ref, "G"),
        ref=record.ref,
        alts=tuple("G"),
    )
    maf_record["COSMIC"] = get_builder("COSMIC", test_scheme, value=None)
    # maf_record["dbSNP_RS"] = get_builder("dbSNP_RS", test_scheme, value="novel")
    found = annotator.annotate(maf_record, vcf_record, var_allele_idx=1)

    # No new dbSNP_RS column without cosmic IDs
    assert found.value == []
    # assert return_nt.dbsnp == None

    ## setup no overlap pos
    vcf_record = get_test_vcf_record(
        chrom=record.chrom,
        pos=101,
        alleles=(record.ref, "G"),
        ref=record.ref,
        alts=tuple("G"),
    )
    maf_record["COSMIC"] = get_builder("COSMIC", test_scheme, value=None)
    # maf_record["dbSNP_RS"] = get_builder("dbSNP_RS", test_scheme, value="novel")
    found = annotator.annotate(maf_record, vcf_record, var_allele_idx=1)

    # No new dbSNP_RS column without cosmic IDs
    assert found.value == []
    # assert return_nt.dbsnp == None


def test_cosmic_del(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_empty_maf_record,
    vcf_gen,
    get_test_vcf_record,
):
    vcf_path = get_test_file("ex2.vcf.gz")
    args = SimpleNamespace(cosmic_vcf=vcf_path)
    annotator = setup_annotator(test_scheme, args)

    gen = vcf_gen("ex2.vcf.gz")

    # deletion
    record = gen.deletion

    ## setup
    vcf_record = get_test_vcf_record(
        chrom=record.chrom,
        pos=record.pos,
        alleles=record.alleles,
        ref=record.ref,
        alts=record.alts,
    )
    maf_record = get_empty_maf_record
    # maf_record["dbSNP_RS"] = get_builder("dbSNP_RS", test_scheme, value="novel")
    found = annotator.annotate(maf_record, vcf_record, var_allele_idx=1)

    assert found.value == ["COSM0001"]
    # assert maf_record["dbSNP_RS"].value == []


def test_cosmic_ins(
    test_scheme,
    setup_annotator,
    get_test_file,
    get_empty_maf_record,
    vcf_gen,
    get_test_vcf_record,
):
    vcf_path = get_test_file("ex2.vcf.gz")
    args = SimpleNamespace(cosmic_vcf=vcf_path)
    annotator = setup_annotator(test_scheme, args)

    gen = vcf_gen("ex2.vcf.gz")

    # insertion
    record = gen.insertion

    ## setup
    vcf_record = get_test_vcf_record(
        chrom=record.chrom,
        pos=record.pos,
        alleles=record.alleles,
        ref=record.ref,
        alts=record.alts,
    )
    maf_record = get_empty_maf_record
    # maf_record["dbSNP_RS"] = get_builder("dbSNP_RS", test_scheme, value="novel")
    found = annotator.annotate(maf_record, vcf_record, var_allele_idx=1)

    assert found.value == ["COSM0002"]
    # assert maf_record["dbSNP_RS"].value == []


# __END__
