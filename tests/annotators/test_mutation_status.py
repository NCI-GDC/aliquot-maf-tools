"""
Tests for the ``aliquotmaf.annotators.MutationStatus`` class.
"""
from collections import OrderedDict

import pytest
from maflib.column_types import MutationStatus
from maflib.column_values import MutationStatusEnum

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


@pytest.fixture
def test_scheme(get_test_scheme):
    coldict = OrderedDict([("Mutation_Status", MutationStatus)])
    return get_test_scheme(coldict)


def test_setup_mutation_status(test_scheme, setup_mutation_status):
    """Test setting up mutation status annotator class"""
    annotator = setup_mutation_status(test_scheme, "MuTect2")
    assert annotator.caller == "MuTect2"


@pytest.mark.parametrize(
    "caller, sample, info, expected",
    [
        ("MuTect2", {"TUMOR": {}}, {}, MutationStatusEnum.Somatic),
        ("GATK4 MuTect2", {"TUMOR": {}}, {}, MutationStatusEnum.Somatic),
        ("SomaticSniper", {"TUMOR": {"SS": 0}}, {}, MutationStatusEnum.NoStatus),
        ("SomaticSniper", {"TUMOR": {"SS": 1}}, {}, MutationStatusEnum.Germline),
        ("SomaticSniper", {"TUMOR": {"SS": 2}}, {}, MutationStatusEnum.Somatic),
        ("SomaticSniper", {"TUMOR": {"SS": 3}}, {}, MutationStatusEnum.LOH),
        ("SomaticSniper", {"TUMOR": {"SS": 4}}, {}, MutationStatusEnum.Unknown),
        ("VarScan2", {}, {"SS": "0"}, MutationStatusEnum.NoStatus),
        ("VarScan2", {}, {"SS": "1"}, MutationStatusEnum.Germline),
        ("VarScan2", {}, {"SS": "2"}, MutationStatusEnum.Somatic),
        ("VarScan2", {}, {"SS": "3"}, MutationStatusEnum.LOH),
        ("VarScan2", {}, {"SS": "5"}, MutationStatusEnum.Unknown),
        ("MuSE", {"TUMOR": {"SS": 0}}, {}, MutationStatusEnum.NoStatus),
        ("MuSE", {"TUMOR": {"SS": 1}}, {}, MutationStatusEnum.Germline),
        ("MuSE", {"TUMOR": {"SS": 2}}, {}, MutationStatusEnum.Somatic),
        ("MuSE", {"TUMOR": {"SS": 3}}, {}, MutationStatusEnum.LOH),
        # ("MuSE", {'TUMOR': {'SS': 4}}, {}, MutationStatusEnum.PostTranscriptional),
        ("MuSE", {"TUMOR": {"SS": 5}}, {}, MutationStatusEnum.Unknown),
    ],
)
def test_mutation_status_annotator(
    test_scheme,
    setup_mutation_status,
    get_test_vcf_record,
    get_empty_maf_record,
    caller,
    sample,
    info,
    expected,
):
    annotator = setup_mutation_status(test_scheme, caller)
    record = get_test_vcf_record(samples=sample, info=info)
    maf_record = annotator.annotate(get_empty_maf_record, record, "TUMOR")
    assert maf_record["Mutation_Status"].value == expected
