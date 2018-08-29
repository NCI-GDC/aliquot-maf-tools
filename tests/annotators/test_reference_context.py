"""
Tests for the ``aliquotmaf.annotators.ReferenceContext`` class.
"""
import pysam
import pytest
from collections import OrderedDict

from maflib.column_types import StringColumn

from aliquotmaf.annotators import ReferenceContext

@pytest.fixture
def setup_annotator():
    created = []
    def _make_annotator(scheme, source, context_size=5):
        curr = ReferenceContext.setup(scheme, source, context_size=context_size)
        created.append(curr)
        return curr

    yield _make_annotator

    for record in created:
        record.shutdown()

@pytest.fixture
def test_scheme(get_test_scheme):
    coldict = OrderedDict([("CONTEXT", StringColumn)])
    return get_test_scheme(coldict)

def test_setup_reference_context(test_scheme, setup_annotator, get_test_file):
    fasta_path = get_test_file('fake_ref.fa')
    annotator = setup_annotator(test_scheme, source=fasta_path) 
    assert annotator.context_size == 5 

def test_reference_context_annotator(test_scheme, setup_annotator, get_test_file,
                                     get_empty_maf_record):
    fasta_path = get_test_file('fake_ref.fa')
    annotator = setup_annotator(test_scheme, source=fasta_path) 

    vcf_path = get_test_file('ex1.vcf.gz')
    vcf_obj = pysam.VariantFile(vcf_path)
    gen = vcf_obj.fetch()

    # simple snp
    record = next(gen)

    ## default context
    maf_record = annotator.annotate(get_empty_maf_record, record) 
    assert maf_record['CONTEXT'].value == 'AGTGGCTCATT'

    ## truncated left context 
    annotator.context_size = 12
    maf_record = annotator.annotate(get_empty_maf_record, record) 
    assert maf_record['CONTEXT'].value == 'CACTAGTGGCTCATTGTAAATG'

    ## truncated both context 
    annotator.context_size = 1575 
    maf_record = annotator.annotate(get_empty_maf_record, record) 
    assert maf_record['CONTEXT'].value == annotator.fa.fetch(region="chr1")

    ## reset
    annotator.context_size = 5

#@pytest.mark.parametrize("caller, sample, info, expected", [
#    ("MuTect2", {'TUMOR': {}}, {}, MutationStatusEnum.Somatic), 
#    ("GATK4 MuTect2", {'TUMOR': {}}, {}, MutationStatusEnum.Somatic), 
#    ("SomaticSniper", {'TUMOR': {'SS': 0}}, {}, MutationStatusEnum.NoStatus),
#    ("SomaticSniper", {'TUMOR': {'SS': 1}}, {}, MutationStatusEnum.Germline),
#    ("SomaticSniper", {'TUMOR': {'SS': 2}}, {}, MutationStatusEnum.Somatic),
#    ("SomaticSniper", {'TUMOR': {'SS': 3}}, {}, MutationStatusEnum.LOH),
#    ("SomaticSniper", {'TUMOR': {'SS': 4}}, {}, MutationStatusEnum.Unknown),
#    ("VarScan2", {}, {'SS': '0'}, MutationStatusEnum.NoStatus),
#    ("VarScan2", {}, {'SS': '1'}, MutationStatusEnum.Germline),
#    ("VarScan2", {}, {'SS': '2'}, MutationStatusEnum.Somatic),
#    ("VarScan2", {}, {'SS': '3'}, MutationStatusEnum.LOH),
#    ("VarScan2", {}, {'SS': '5'}, MutationStatusEnum.Unknown),
#    ("MuSE", {'TUMOR': {'SS': 0}}, {}, MutationStatusEnum.NoStatus),
#    ("MuSE", {'TUMOR': {'SS': 1}}, {}, MutationStatusEnum.Germline),
#    ("MuSE", {'TUMOR': {'SS': 2}}, {}, MutationStatusEnum.Somatic),
#    ("MuSE", {'TUMOR': {'SS': 3}}, {}, MutationStatusEnum.LOH),
#    #("MuSE", {'TUMOR': {'SS': 4}}, {}, MutationStatusEnum.PostTranscriptional),
#    ("MuSE", {'TUMOR': {'SS': 5}}, {}, MutationStatusEnum.Unknown),
#])
#def test_mutation_status_annotator(test_scheme, setup_mutation_status, 
#                                   get_test_vcf_record, get_empty_maf_record, 
#                                   caller, sample, info, expected):
#    annotator = setup_mutation_status(test_scheme, caller)
#    record = get_test_vcf_record(samples=sample, info=info)
#    maf_record = annotator.annotate(get_empty_maf_record, record, 'TUMOR')
#    assert maf_record['Mutation_Status'].value == expected 
