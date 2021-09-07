import logging
import os
from collections import namedtuple

import pysam
import pytest
from maflib.record import MafRecord
from maflib.schemes import MafScheme
from maflib.validation import ValidationStringency

log = logging.getLogger()
log.setLevel(logging.INFO)

GenericVcfRecord = namedtuple(
    "GenericVcfRecord",
    ["chrom", "pos", "stop", "alleles", "ref", "alts", "samples", "info"],
)

GenericVcfObject = namedtuple("GenericVcfObject", ["header"])


@pytest.fixture
def get_test_scheme():
    """
    Generates a very basice scheme used for testing.
    """

    def generate_scheme(
        coldict,
        name="_TestScheme",
        version="test-version",
        annotation="test-annotation",
    ):
        tpe = type(str(name), (MafScheme,), {})
        setattr(tpe, "version", classmethod(lambda cls: version))
        setattr(tpe, "annotation_spec", classmethod(lambda cls: annotation))
        setattr(tpe, "__column_dict__", classmethod(lambda cls: coldict))
        return tpe()

    return generate_scheme


@pytest.fixture
def get_test_vcf_record():
    """
    Fixture to mock a vcf record.
    """

    def generate_vcf_record(
        chrom=None,
        pos=None,
        stop=None,
        alleles=(),
        ref=None,
        alts=(),
        samples={},
        info={},
    ):
        tpe = GenericVcfRecord(
            chrom=chrom,
            pos=pos,
            stop=stop,
            alleles=alleles,
            ref=ref,
            alts=alts,
            samples=samples,
            info=info,
        )

        return tpe

    return generate_vcf_record


@pytest.fixture
def get_empty_maf_record():
    """
    Generates an empty MAF record.
    """

    def _get_empty_maf_record(line_number=0, stringency=ValidationStringency.Strict):
        return MafRecord(line_number=line_number, validation_stringency=stringency)

    return _get_empty_maf_record()


@pytest.fixture
def get_test_file():
    """
    Loads a file from data.
    """

    def _get_test_file(name):
        curr = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(curr, "data/{0}".format(name))

    return _get_test_file


@pytest.fixture
def get_test_vcf_header():
    """
    Generate VCF header.
    """

    def _get_test_vcf_header(meta=None, samples=None):
        hdr = pysam.VariantHeader()
        if meta:
            if isinstance(meta, list):
                for rec in meta:
                    hdr.add_meta(**rec)
            else:
                hdr.add_meta(**meta)
                # hdr.add_meta(key=meta['key'], value=meta.get('value'), items=meta.get('items'))

        if samples:
            if isinstance(samples, list):
                for sample in samples:
                    hdr.add_sample(sample)
            else:
                hdr.add_sample(samples)

        res = GenericVcfObject(header=hdr)
        return res

    return _get_test_vcf_header
