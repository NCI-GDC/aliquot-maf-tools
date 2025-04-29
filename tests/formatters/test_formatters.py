"""
Tests for functions in `aliquotmaf.converters.formatters` which are used to
format particular MAF column values.
"""
import pytest

import aliquotmaf.converters.formatters as Formatters
import aliquotmaf.subcommands.vcf_to_aliquot.extractors as Extractors
from aliquotmaf.constants import variant_callers


@pytest.mark.parametrize(
    "genotype, ref_allele, var_allele, position, alleles, caller_id, expected",
    [
        ({"GT": {}}, "A", "T", 100, ("A", "T"), variant_callers.MUTECT2, ("A", "T")),
        (
            {"GT": (0, 1), "AD": (1, 9), "DP": 10},
            "A",
            "T",
            100,
            ("A", "T"),
            variant_callers.MUTECT2,
            ("A", "T"),
        ),
        (
            {"GT": (0, 0), "AD": (1, 9), "DP": 10},
            "A",
            "T",
            100,
            ("A", "T"),
            variant_callers.MUTECT2,
            ("A", "A"),
        ),
        (
            {"GT": (0, 1), "AD": (1, 9), "DP": 10},
            "AT",
            "A",
            100,
            ("AT", "A"),
            variant_callers.MUTECT2,
            ("T", "-"),
        ),
    ],
)
def test_format_alleles(
    genotype, ref_allele, var_allele, position, alleles, caller_id, expected
):
    """
    Tests the function which formats alleles for MAF records
    """
    idx = Extractors.VariantAlleleIndexExtractor.extract(genotype)
    gt, depths = Extractors.GenotypeAndDepthsExtractor.extract(
        idx, genotype, alleles, caller_id
    )
    loc = Extractors.LocationDataExtractor.extract(
        ref_allele, var_allele, position, alleles
    )

    a1, a2 = Formatters.format_alleles(
        genotype=gt,
        alleles=loc["alleles"],
        defaults=[loc["ref_allele"], loc["var_allele"]],
    )

    assert (a1, a2) == expected


# @pytest.mark.parametrize(
#     "genotype, alleles, expected",
#     [
#         ({"GT": {}}, ("A", "T"), (0, None, None)),
#         ({"GT": (0, 1), "AD": (1, 9), "DP": 10}, ("A", "T"), (10, 1, 9)),
#         ({"GT": (0, 2), "AD": (1, 0, 9), "DP": 10}, ("A", "G", "T"), (10, 1, 9)),
#     ],
# )
# def test_format_depths(genotype, alleles, expected):
#     """
#     Tests the function which formats depths for the MAF.
#     """
#     idx = Extractors.VariantAlleleIndexExtractor.extract(genotype)
#     gt, depths = Extractors.GenotypeAndDepthsExtractor.extract(idx, genotype, alleles)

#     dp, ref_ct, alt_ct = Formatters.format_depths(
#         genotype=gt, depths=depths, var_allele_idx=idx, default_total_dp=0
#     )
#     assert (dp, ref_ct, alt_ct) == expected
