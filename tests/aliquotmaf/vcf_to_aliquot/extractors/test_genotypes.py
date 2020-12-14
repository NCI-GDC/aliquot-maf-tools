"""
Tests the extractors in aliquotmaf.subcommands.vcf_to_protected.extractors.genotypes
"""
import pytest

from aliquotmaf.vcf_to_aliquot.extractors.genotypes import (
    GenotypeAndDepthsExtractor,
    VariantAlleleIndexExtractor,
)


@pytest.mark.parametrize(
    "genotype, expected",
    [
        ({"GT": (0, 0)}, 1),
        ({"GT": (0, 1)}, 1),
        ({"GT": (1, 0)}, 1),
        ({"GT": (2, 1)}, 2),
        ({"GT": (1, 2)}, 1),
        ({"GT": (None, None)}, 1),
    ],
)
def test_variant_allele_index_extractor(genotype, expected):
    """
    Tests the extraction of the variant allele index from the vcf
    """
    idx = VariantAlleleIndexExtractor.extract(genotype)
    assert idx == expected


@pytest.mark.parametrize(
    "genotype, alleles, expected_gt, expected_dp",
    [
        ({"GT": None}, (), {}, []),
        (
            {"GT": (0, 0), "DP": 10, "AD": (9, 1)},
            ("A", "T"),
            {"GT": (0, 0), "DP": 10, "AD": (9, 1)},
            [9, 1],
        ),
        (
            {"GT": (0, 0), "DP": 10, "AD": 8, "RD": 2},
            ("A", "T"),
            {"GT": (0, 0), "DP": 10, "AD": (2, 8)},
            [2, 8],
        ),
        (
            {"GT": (0, 0), "DP": 10, "AD": (8), "RD": 2},
            ("A", "T"),
            {"GT": (0, 0), "DP": 10, "AD": (2, 8)},
            [2, 8],
        ),
        (
            {"GT": (0, 0), "DP": 10, "BCOUNT": (8, 0, 0, 2)},
            ("A", "T"),
            {"GT": (0, 0), "DP": 10, "AD": (8, 2)},
            [8, 2],
        ),
        (
            {"GT": (0, 0), "DP": 10, "AD": (11, 1)},
            ("A", "T"),
            {"GT": (0, 0), "DP": 12, "AD": (11, 1)},
            [11, 1],
        ),
        (
            {"GT": (0, 0), "DP": 10, "AD": (1, 11)},
            ("A", "T"),
            {"GT": (0, 0), "DP": 12, "AD": (1, 11)},
            [1, 11],
        ),
        (
            {"GT": (0, 0), "DP": 10, "AD": (1, 10)},
            ("A", "T"),
            {"GT": (0, 0), "DP": 11, "AD": (1, 10)},
            [1, 10],
        ),
    ],
)
def test_genotype_and_depths_extractor(genotype, alleles, expected_gt, expected_dp):
    """
    Tests the extraction of the genotype and depths data from the vcf
    """
    idx = VariantAlleleIndexExtractor.extract(genotype)
    new_gt, depths = GenotypeAndDepthsExtractor.extract(idx, genotype, alleles)
    assert new_gt == expected_gt
    assert depths == expected_dp
