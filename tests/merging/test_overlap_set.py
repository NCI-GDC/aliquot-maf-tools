"""
Tests associated with the OverlapSet class.
"""


def test_overlap_set_basic_singleton(test_input_scheme, overlapped_records_generate):
    """
    Test simple parsing of overlaps with singleton
    """
    callers = ["MuTect2"]
    maf_lines = [
        "chr1\t1\t1\tSNP\tA\tC\tA\tC\tA\tA\t10\t2\t8\t8\t8\t0\t\t\n",
    ]

    record = overlapped_records_generate(test_input_scheme, maf_lines, callers)

    assert record.is_singleton()

    assert ["MuTect2"] == record.callers

    assert ("SNP",) == record.variant_types

    assert "1:1:C" in record.locus_allele_map
    assert len(record.locus_allele_map) == 1
    assert len(record.locus_allele_map["1:1:C"]) == 1

    assert ("MuTect2", "SNP") in record.caller_type_map
    assert len(record.caller_type_map) == 1

    assert record.all_single_record() is True


def test_overlap_set_basic_a(test_input_scheme, overlapped_records_generate):
    """
    Test simple parsing of overlaps
    """
    callers = ["MuTect2", "MuSE"]
    maf_lines = [
        "chr1\t1\t1\tSNP\tA\tC\tA\tC\tA\tA\t10\t2\t8\t8\t8\t0\t\t\n",
        "chr1\t1\t1\tSNP\tA\tC\tA\tC\tA\tA\t20\t2\t18\t8\t8\t0\t\t\n",
    ]

    record = overlapped_records_generate(test_input_scheme, maf_lines, callers)

    assert not record.is_singleton()

    assert ["MuSE", "MuTect2"] == record.callers

    assert ("SNP",) == record.variant_types

    assert "1:1:C" in record.locus_allele_map
    assert len(record.locus_allele_map) == 1
    assert len(record.locus_allele_map["1:1:C"]) == 2

    assert ("MuTect2", "SNP") in record.caller_type_map and (
        "MuSE",
        "SNP",
    ) in record.caller_type_map
    assert len(record.caller_type_map) == 2

    assert record.all_single_record() is True


def test_overlap_set_basic_b(test_input_scheme, overlapped_records_generate):
    """
    Test simple parsing of overlaps when one caller has no overlaps
    """
    callers = ["MuTect2", "MuSE", "SomaticSniper"]
    maf_lines = [
        "chr1\t1\t1\tSNP\tA\tC\tA\tC\tA\tA\t10\t2\t8\t8\t8\t0\t\t\n",
        "chr1\t1\t1\tSNP\tA\tC\tA\tC\tA\tA\t20\t2\t18\t8\t8\t0\t\t\n",
        "",
    ]

    record = overlapped_records_generate(test_input_scheme, maf_lines, callers)

    assert not record.is_singleton()

    assert ["MuSE", "MuTect2"] == record.callers

    assert ("SNP",) == record.variant_types

    assert "1:1:C" in record.locus_allele_map
    assert len(record.locus_allele_map) == 1
    assert len(record.locus_allele_map["1:1:C"]) == 2

    assert ("MuTect2", "SNP") in record.caller_type_map and (
        "MuSE",
        "SNP",
    ) in record.caller_type_map
    assert len(record.caller_type_map) == 2

    assert record.all_single_record() is True


def test_overlap_set_basic_c(test_input_scheme, overlapped_records_generate):
    """
    Test simple parsing of overlaps with 3 callers and 2 unique alleles
    """
    callers = ["MuTect2", "MuSE", "SomaticSniper"]
    maf_lines = [
        "chr1\t1\t1\tSNP\tA\tC\tA\tC\tA\tA\t10\t2\t8\t8\t8\t0\t\t\n",
        "chr1\t1\t1\tSNP\tA\tC\tA\tC\tA\tA\t20\t2\t18\t8\t8\t0\t\t\n",
        "chr1\t1\t1\tSNP\tA\tT\tA\tT\tA\tA\t20\t2\t18\t8\t8\t0\t\t\n",
    ]

    record = overlapped_records_generate(test_input_scheme, maf_lines, callers)

    assert not record.is_singleton()

    assert ["MuSE", "MuTect2", "SomaticSniper"] == record.callers

    assert ("SNP",) == record.variant_types

    assert "1:1:C" in record.locus_allele_map
    assert "1:1:T" in record.locus_allele_map
    assert len(record.locus_allele_map) == 2
    assert len(record.locus_allele_map["1:1:C"]) == 2
    assert len(record.locus_allele_map["1:1:T"]) == 1

    assert (
        ("MuTect2", "SNP") in record.caller_type_map
        and ("MuSE", "SNP") in record.caller_type_map
        and ("SomaticSniper", "SNP") in record.caller_type_map
    )
    assert len(record.caller_type_map) == 3

    assert record.all_single_record() is True


def test_overlap_set_basic_d(test_input_scheme, overlapped_records_generate):
    """
    Test simple parsing of overlaps with 3 callers, 1 caller with multiple
    overlaps between SNP and INS.
    """
    callers = ["MuTect2", "MuSE", "SomaticSniper"]
    maf_lines = [
        [
            "chr1\t1\t1\tSNP\tA\tC\tA\tC\tA\tA\t10\t2\t8\t8\t8\t0\t\t\n",
            "chr1\t1\t1\tINS\t-\tCC\t-\tCC\t-\t-\t10\t2\t8\t8\t8\t0\t\t\n",
        ],
        "chr1\t1\t1\tSNP\tA\tC\tA\tC\tA\tA\t20\t2\t18\t8\t8\t0\t\t\n",
        "chr1\t1\t1\tSNP\tA\tT\tA\tT\tA\tA\t20\t2\t18\t8\t8\t0\t\t\n",
    ]

    record = overlapped_records_generate(test_input_scheme, maf_lines, callers)

    assert not record.is_singleton()

    assert ["MuSE", "MuTect2", "SomaticSniper"] == record.callers

    assert ("INS", "SNP") == record.variant_types

    assert "1:1:C" in record.locus_allele_map
    assert "1:1:CC" in record.locus_allele_map
    assert "1:1:T" in record.locus_allele_map
    assert len(record.locus_allele_map) == 3
    assert len(record.locus_allele_map["1:1:C"]) == 2
    assert len(record.locus_allele_map["1:1:T"]) == 1
    assert len(record.locus_allele_map["1:1:CC"]) == 1

    assert (
        ("MuTect2", "SNP") in record.caller_type_map
        and ("MuSE", "SNP") in record.caller_type_map
        and ("SomaticSniper", "SNP") in record.caller_type_map
        and ("MuTect2", "INS") in record.caller_type_map
    )

    assert len(record.caller_type_map) == 4

    assert record.all_single_record() is False
