"""
Tests of the aliquotmaf.merging.record_merger.impl.v1_0.MafRecordMerger_1_0_0 class
"""

from aliquotmaf.merging.record_merger.impl.v1_0 import MafRecordMerger_1_0_0

def test_record_merge_1(test_input_scheme, test_output_scheme, overlapped_records_generate):
    callers = ['mutect2']
    maf_lines = [
        'chr1\t1\t1\tSNP\tA\tC\tA\tC\tA\tA\t10\t2\t8\t8\t8\t0\t\t\n',
    ]

    record = overlapped_records_generate(
        test_input_scheme,
        maf_lines,
        callers)

    merger = MafRecordMerger_1_0_0(test_output_scheme)

    results = merger.merge_records(record)
    assert len(results) == 1
    result = results[0]

    assert result["callers"].value == ['mutect2']

def test_record_merge_2(test_input_scheme, test_output_scheme, overlapped_records_generate):
    callers = ['mutect2']
    maf_lines = [
        [
            'chr1\t1\t1\tSNP\tA\tC\tA\tC\tA\tA\t10\t2\t8\t8\t8\t0\t\t\n',
            'chr1\t1\t1\tINS\t-\tCC\t-\tCC\t-\t-\t10\t2\t8\t8\t8\t0\t\t\n',
        ]
    ]

    record = overlapped_records_generate(
        test_input_scheme,
        maf_lines,
        callers)

    merger = MafRecordMerger_1_0_0(test_output_scheme)

    results = merger.merge_records(record)
    assert len(results) == 2

    for result in results:
        assert result["callers"].value == ['mutect2']

def test_record_merge_3(test_input_scheme, test_output_scheme, overlapped_records_generate):
    callers = ['mutect2', 'muse']
    maf_lines = [
        'chr1\t1\t1\tSNP\tA\tC\tA\tC\tA\tA\t10\t2\t8\t8\t8\t0\tA\tgdc_pon\n',
        'chr1\t1\t1\tSNP\tA\tC\tC\tC\tA\tA\t12\t2\t10\t20\t20\t0\tB\tgdc_pon;off_target\n',
    ]

    record = overlapped_records_generate(
        test_input_scheme,
        maf_lines,
        callers)

    merger = MafRecordMerger_1_0_0(test_output_scheme)

    results = merger.merge_records(record)
    assert len(results) == 1

    for result in results:
        assert result['callers'].value == ['muse', 'mutect2']
        assert result['extra'].value == 'A'
        assert result['GDC_FILTER'].value == ['gdc_pon', 'off_target']
        assert result['Tumor_Seq_Allele1'].value == 'A'
        assert result['Tumor_Seq_Allele2'].value == 'C'
        assert result['t_depth'].value == 11
        assert result['t_ref_count'].value == 2
        assert result['t_alt_count'].value == 9
        assert result['n_depth'].value == 14
        assert result['n_ref_count'].value == 14
        assert result['n_alt_count'].value == 0

def test_record_merge_4(test_input_scheme, test_output_scheme, overlapped_records_generate):
    callers = ['mutect2', 'muse', 'pindel']
    maf_lines = [
        'chr1\t1\t1\tSNP\tA\tC\tA\tC\tA\tA\t10\t2\t8\t8\t8\t0\tA\t\n',
        'chr1\t1\t1\tSNP\tA\tC\tC\tC\tA\tA\t12\t2\t10\t20\t20\t0\tB\tgdc_pon\n',
        'chr1\t1\t2\tDNP\tAT\tCG\tAT\tCG\tAT\tAT\t8\t2\t6\t10\t10\t0\tC\toff_target\n'
    ]

    record = overlapped_records_generate(
        test_input_scheme,
        maf_lines,
        callers)

    merger = MafRecordMerger_1_0_0(test_output_scheme)

    results = merger.merge_records(record)
    assert len(results) == 1
    for result in results:
        assert result['callers'].value == ['muse*', 'mutect2*', 'pindel']
        assert result['extra'].value == 'C'
        assert result['GDC_FILTER'].value == ['off_target']
        assert result['Variant_Type'].value.value == 'DNP'
        assert result['Reference_Allele'].value == 'AT'
        assert result['Tumor_Seq_Allele1'].value == 'AT'
        assert result['Tumor_Seq_Allele2'].value == 'CG'
        assert result['t_depth'].value == 8 
        assert result['t_ref_count'].value == 2
        assert result['t_alt_count'].value == 6
        assert result['n_depth'].value == 10 
        assert result['n_ref_count'].value == 10
        assert result['n_alt_count'].value == 0
