import pytest
from collections import OrderedDict

from maflib.column_types import MutationStatus
from aliquotmaf.annotators import MutationStatus as MutationStatusAnnotator

@pytest.fixture
def setup_mutation_status():
    created = []
    def _make_mutation_status(scheme, caller):
        curr = MutationStatusAnnotator.setup(scheme, caller)
        created.append(curr)
        return curr

    yield _make_mutation_status

    for record in created:
        record.shutdown()

#from maflib.validation import ValidationStringency
#from maflib.record import MafRecord

@pytest.fixture
def test_scheme(get_test_scheme):
    coldict = OrderedDict([("Mutation_Status", MutationStatus)])
    return get_test_scheme(coldict)

#@pytest.fixture
#def get_empty_maf_record():
#
#    def _get_empty_maf_record(line_number=0, stringency=ValidationStringency.Strict):
#        return MafRecord(line_number=line_number, validation_stringency=stringency)
#
#    return _get_empty_maf_record()
