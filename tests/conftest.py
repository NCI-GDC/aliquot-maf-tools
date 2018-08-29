import logging
import pytest
from collections import OrderedDict, namedtuple

from maflib.schemes import MafScheme
from maflib.validation import ValidationStringency
from maflib.record import MafRecord

log = logging.getLogger()
log.setLevel(logging.INFO)

GenericVcfRecord = namedtuple("GenericVcfRecord", ["samples", "info"])

@pytest.fixture
def get_test_scheme():
    def generate_scheme(coldict, name="_TestScheme", version="test-version", annotation="test-annotation"):
        tpe = type(str(name), (MafScheme,), {})
        setattr(tpe, "version", classmethod(lambda cls: version))
        setattr(tpe, "annotation_spec", classmethod(lambda cls: annotation))
        setattr(tpe, "__column_dict__", classmethod(lambda cls: coldict))
        return tpe() 

    return generate_scheme

@pytest.fixture
def get_test_vcf_record():
    def generate_vcf_record(samples={}, info={}):
        tpe = GenericVcfRecord(
            samples=samples,
            info=info
        )

        return tpe

    return generate_vcf_record

@pytest.fixture
def get_empty_maf_record():

    def _get_empty_maf_record(line_number=0, stringency=ValidationStringency.Strict):
        return MafRecord(line_number=line_number, validation_stringency=stringency)

    return _get_empty_maf_record()
