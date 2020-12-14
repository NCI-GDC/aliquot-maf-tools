"""
Tests the extractors in aliquotmaf.subcommands.vcf_to_aliquot.extractors.location
"""
import pytest

from aliquotmaf.vcf_to_aliquot.extractors.location import LocationDataExtractor

# VariantAlleleIndexExtractor -> GenotypeAndDepthsExtractor -> LocationDataExtractor -> EffectsExtractor -> SelectOneEffectExtractor -> PopulationFrequencyExtractor -> VariantClassExtractor


@pytest.mark.parametrize(
    "ref_allele, var_allele, position, alleles, expected",
    [
        (
            "A",
            "T",
            100,
            ("A", "T"),
            {
                "alleles": ("A", "T"),
                "ref_allele": "A",
                "var_allele": "T",
                "start": 100,
                "stop": 100,
                "var_type": "SNP",
                "inframe": None,
            },
        ),
        (
            "AC",
            "TA",
            100,
            ("AC", "TA"),
            {
                "alleles": ("AC", "TA"),
                "ref_allele": "AC",
                "var_allele": "TA",
                "start": 100,
                "stop": 101,
                "var_type": "DNP",
                "inframe": None,
            },
        ),
        (
            "ACT",
            "TAC",
            100,
            ("ACT", "TAC"),
            {
                "alleles": ("ACT", "TAC"),
                "ref_allele": "ACT",
                "var_allele": "TAC",
                "start": 100,
                "stop": 102,
                "var_type": "TNP",
                "inframe": None,
            },
        ),
        (
            "ACTT",
            "TACC",
            100,
            ("ACTT", "TACC"),
            {
                "alleles": ("ACTT", "TACC"),
                "ref_allele": "ACTT",
                "var_allele": "TACC",
                "start": 100,
                "stop": 103,
                "var_type": "ONP",
                "inframe": None,
            },
        ),
    ],
)
def test_location_data_extractor_nps(
    ref_allele, var_allele, position, alleles, expected
):
    """
    Tests the normalization of alleles and postions into the MAF format and
    the detection of variant types for SNP/MNPs
    """
    res = LocationDataExtractor.extract(ref_allele, var_allele, position, alleles)
    assert res == expected


@pytest.mark.parametrize(
    "ref_allele, var_allele, position, alleles, expected",
    [
        (
            "A",
            "AT",
            100,
            ("A", "AT"),
            {
                "alleles": ["-", "T"],
                "ref_allele": "-",
                "var_allele": "T",
                "start": 100,
                "stop": 101,
                "var_type": "INS",
                "inframe": False,
            },
        ),
        (
            "A",
            "ATATATT",
            100,
            ("A", "ATATATT"),
            {
                "alleles": ["-", "TATATT"],
                "ref_allele": "-",
                "var_allele": "TATATT",
                "start": 100,
                "stop": 101,
                "var_type": "INS",
                "inframe": True,
            },
        ),
    ],
)
def test_location_data_extractor_ins(
    ref_allele, var_allele, position, alleles, expected
):
    """
    Tests the normalization of alleles and postions into the MAF format and
    the detection of variant types for insertions
    """
    res = LocationDataExtractor.extract(ref_allele, var_allele, position, alleles)
    assert res == expected


@pytest.mark.parametrize(
    "ref_allele, var_allele, position, alleles, expected",
    [
        (
            "AT",
            "A",
            100,
            ("AT", "A"),
            {
                "alleles": ["T", "-"],
                "ref_allele": "T",
                "var_allele": "-",
                "start": 101,
                "stop": 101,
                "var_type": "DEL",
                "inframe": False,
            },
        ),
        (
            "ATT",
            "A",
            100,
            ("ATT", "A"),
            {
                "alleles": ["TT", "-"],
                "ref_allele": "TT",
                "var_allele": "-",
                "start": 101,
                "stop": 102,
                "var_type": "DEL",
                "inframe": False,
            },
        ),
        (
            "ATTT",
            "A",
            100,
            ("ATTT", "A"),
            {
                "alleles": ["TTT", "-"],
                "ref_allele": "TTT",
                "var_allele": "-",
                "start": 101,
                "stop": 103,
                "var_type": "DEL",
                "inframe": True,
            },
        ),
    ],
)
def test_location_data_extractor_del(
    ref_allele, var_allele, position, alleles, expected
):
    """
    Tests the normalization of alleles and postions into the MAF format and
    the detection of variant types for deletions
    """
    res = LocationDataExtractor.extract(ref_allele, var_allele, position, alleles)
    assert res == expected


@pytest.mark.parametrize(
    "ref_allele, var_allele, position, alleles, expected",
    [
        (
            "AT",
            "ACC",
            100,
            ("AT", "ACC"),
            {
                "alleles": ["T", "CC"],
                "ref_allele": "T",
                "var_allele": "CC",
                "start": 101,
                "stop": 101,
                "var_type": "INS",
                "inframe": False,
            },
        ),
        (
            "ATT",
            "AC",
            100,
            ("ATT", "AC"),
            {
                "alleles": ["TT", "C"],
                "ref_allele": "TT",
                "var_allele": "C",
                "start": 101,
                "stop": 102,
                "var_type": "DEL",
                "inframe": False,
            },
        ),
    ],
)
def test_location_data_extractor_indel(
    ref_allele, var_allele, position, alleles, expected
):
    """
    Tests the normalization of alleles and postions into the MAF format and
    the detection of variant types for indels.
    """
    res = LocationDataExtractor.extract(ref_allele, var_allele, position, alleles)
    assert res == expected
