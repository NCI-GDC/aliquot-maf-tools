"""
Tests for the ``aliquotmaf.annotators.MutationStatus`` class.
"""
import pytest

from maflib.column_values import MutationStatusEnum
from aliquotmaf.annotators import MutationStatus

@pytest.fixture
def setup_mutation_status():
    created = []
    def _make_mutation_status(scheme, caller):
        curr = MutationStatus.setup(scheme, caller) 
        created.append(curr)
        return curr

    yield _make_mutation_status

    for record in created:
        record.shutdown()

def test_setup_mutation_status(get_test_scheme, setup_mutation_status):
    """Test setting up mutation status annotator class"""
    annotator = setup_mutation_status(get_test_scheme, 'MuTect2')
    assert annotator.caller == 'MuTect2'

def test_mutation_status_mutect2(get_test_scheme, setup_mutation_status, get_test_vcf_record, get_empty_maf_record):
    annotator = setup_mutation_status(get_test_scheme, 'MuTect2')
    maf_record = annotator.annotate(get_empty_maf_record, get_test_vcf_record, 'TUMOR')
    assert maf_record['Mutation_Status'].value == MutationStatusEnum.Somatic 
