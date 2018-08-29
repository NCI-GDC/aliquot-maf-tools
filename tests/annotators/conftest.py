import pytest
from collections import OrderedDict

from maflib.schemes import MafScheme
from maflib.validation import ValidationStringency
from maflib.column_types import StringColumn, FloatColumn, MutationStatus
from maflib.record import MafRecord

class _TestScheme(MafScheme):
    @classmethod
    def version(cls):
        return "test-version"

    @classmethod
    def annotation_spec(cls):
        return "test-annotation"

    @classmethod
    def __column_dict__(cls):
        return OrderedDict([("Mutation_Status", MutationStatus),
                            ("float", FloatColumn),
                            ("str2", StringColumn)])

class _TestVcfRecord:
    def __init__(self):
        self.samples = {'TUMOR': {}} 
        self.info = {}

@pytest.fixture
def get_test_scheme():
    def generate_scheme():
        return _TestScheme()
    return generate_scheme()

@pytest.fixture
def get_test_vcf_record():
    def generate_vcf_record():
        return _TestVcfRecord()

    return generate_vcf_record()

@pytest.fixture
def get_empty_maf_record():

    def _get_empty_maf_record(line_number=0, stringency=ValidationStringency.Strict):
        return MafRecord(line_number=line_number, validation_stringency=stringency)

    return _get_empty_maf_record()
