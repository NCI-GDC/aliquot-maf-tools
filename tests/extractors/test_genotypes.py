"""
Tests the extractors in aliquotmaf.subcommands.vcf_to_protected.extractors.genotypes
"""
from unittest import TestCase
from unittest.mock import patch

import pytest

from aliquotmaf.constants import variant_callers
from aliquotmaf.subcommands.vcf_to_aliquot.extractors.genotypes import (
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


class TestGenotypeExtractor(TestCase):
    @patch(
        "aliquotmaf.subcommands.vcf_to_aliquot.extractors.genotypes.GenotypeAndDepthsExtractor._dispatch_extractor",
        autospec=True,
        return_value=([5, 6], 11),
    )
    def test_extract_empty_gt(self, dispatch_extractor):
        """
        Test case where GT is present
        """
        # Prepare
        idx = 1
        genotype = {'GT': ''}
        alleles = []
        caller_id = 'somecaller'
        # Test
        new_gt, depths = GenotypeAndDepthsExtractor.extract(
            idx, genotype, alleles, caller_id
        )
        # Expected results
        expected_newgt = {}
        expected_depths = []
        # Validate results
        dispatch_extractor.assert_not_called()
        assert new_gt == expected_newgt
        assert depths == expected_depths

    @patch(
        "aliquotmaf.subcommands.vcf_to_aliquot.extractors.genotypes.GenotypeAndDepthsExtractor._dispatch_extractor",
        autospec=True,
        return_value=([5, 6], 11),
    )
    def test_extract_with_gt(self, dispatch_extractor):
        """
        Test case where GT is present
        """
        # Prepare
        idx = 1
        genotype = {'GT': '0/1'}
        alleles = []
        caller_id = 'somecaller'
        # Test
        new_gt, depths = GenotypeAndDepthsExtractor.extract(
            idx, genotype, alleles, caller_id
        )
        # Expected results
        expected_newgt = {
            'AD': (5, 6),
            'GT': '0/1',
            'DP': 11,
        }
        expected_depths = [5, 6]
        # Validate results
        dispatch_extractor.assert_called_once_with(idx, genotype, alleles, caller_id)
        assert new_gt == expected_newgt
        assert depths == expected_depths

    @patch(
        "aliquotmaf.subcommands.vcf_to_aliquot.extractors.genotypes.GenotypeAndDepthsExtractor._extract_caveman",
        autospec=True,
        return_value=([5, 6], 11),
    )
    def test_extract_caveman(self, extract_fn):
        """
        Test dispatch CAVEMAN
        """
        # Prepare
        caller = variant_callers.CAVEMAN
        # Test
        newgt, depths = GenotypeAndDepthsExtractor._dispatch_extractor(
            1, {'GT': '0/1'}, [], caller
        )
        # Validate results
        extract_fn.assert_called_once_with(1, {'GT': '0/1'}, [])

    @patch(
        "aliquotmaf.subcommands.vcf_to_aliquot.extractors.genotypes.GenotypeAndDepthsExtractor._extract_mutect2",
        autospec=True,
        return_value=([5, 6], 11),
    )
    def test_extract_gatk4_mutect2(self, extract_fn):
        """
        Test dispatch GATK4 Mutect2
        """
        # Prepare
        caller = variant_callers.GATK4_MUTECT2
        # Test
        newgt, depths = GenotypeAndDepthsExtractor._dispatch_extractor(
            1, {'GT': '0/1'}, [], caller
        )
        # Validate results
        extract_fn.assert_called_once_with({'GT': '0/1'})

    @patch(
        "aliquotmaf.subcommands.vcf_to_aliquot.extractors.genotypes.GenotypeAndDepthsExtractor._extract_mutect2",
        autospec=True,
        return_value=([5, 6], 11),
    )
    def test_extract_mutect2(self, extract_fn):
        """
        Test dispatch Mutect2
        """
        # Prepare
        caller = variant_callers.MUTECT2
        # Test
        newgt, depths = GenotypeAndDepthsExtractor._dispatch_extractor(
            1, {'GT': '0/1'}, [], caller
        )
        # Validate results
        extract_fn.assert_called_once_with({'GT': '0/1'})

    @patch(
        "aliquotmaf.subcommands.vcf_to_aliquot.extractors.genotypes.GenotypeAndDepthsExtractor._extract_mutect2",
        autospec=True,
        return_value=([5, 6], 11),
    )
    def test_extract_gatk4_mutect2_pair(self, extract_fn):
        """
        Test dispatch GATK4 Mutect2 Pair
        """
        # Prepare
        caller = variant_callers.GATK4_MUTECT2_PAIR
        # Test
        newgt, depths = GenotypeAndDepthsExtractor._dispatch_extractor(
            1, {'GT': '0/1'}, [], caller
        )
        # Validate results
        extract_fn.assert_called_once_with({'GT': '0/1'})

    @patch(
        "aliquotmaf.subcommands.vcf_to_aliquot.extractors.genotypes.GenotypeAndDepthsExtractor._extract_somaticsniper",
        autospec=True,
        return_value=([5, 6], 11),
    )
    def test_extract_somtaticsniper(self, extract_fn):
        """
        Test dispatch Somatic Sniper
        """
        # Prepare
        caller = variant_callers.SOMATIC_SNIPER
        # Test
        newgt, depths = GenotypeAndDepthsExtractor._dispatch_extractor(
            1, {'GT': '0/1'}, [], caller
        )
        # Validate results
        extract_fn.assert_called_once_with(1, {'GT': '0/1'}, [])

    @patch(
        "aliquotmaf.subcommands.vcf_to_aliquot.extractors.genotypes.GenotypeAndDepthsExtractor._extract_muse",
        autospec=True,
        return_value=([5, 6], 11),
    )
    def test_extract_muse(self, extract_fn):
        """
        Test dispatch Muse
        """
        # Prepare
        caller = variant_callers.MUSE
        # Test
        newgt, depths = GenotypeAndDepthsExtractor._dispatch_extractor(
            1, {'GT': '0/1'}, [], caller
        )
        # Validate results
        extract_fn.assert_called_once_with({'GT': '0/1'})

    @patch(
        "aliquotmaf.subcommands.vcf_to_aliquot.extractors.genotypes.GenotypeAndDepthsExtractor._extract_varscan2",
        autospec=True,
        return_value=([5, 6], 11),
    )
    def test_extract_varscan2(self, extract_fn):
        """
        Test dispatch Varscan 2
        """
        # Prepare
        caller = variant_callers.VARSCAN2
        # Test
        newgt, depths = GenotypeAndDepthsExtractor._dispatch_extractor(
            1, {'GT': '0/1'}, [], caller
        )
        # Validate results
        extract_fn.assert_called_once_with(1, {'GT': '0/1'}, [])

    @patch(
        "aliquotmaf.subcommands.vcf_to_aliquot.extractors.genotypes.GenotypeAndDepthsExtractor._extract_pindel",
        autospec=True,
        return_value=([5, 6], 11),
    )
    def test_extract_pindel(self, extract_fn):
        """
        Test dispatch Pindel
        """
        # Prepare
        caller = variant_callers.PINDEL
        # Test
        newgt, depths = GenotypeAndDepthsExtractor._dispatch_extractor(
            1, {'GT': '0/1'}, [], caller
        )
        # Validate results
        extract_fn.assert_called_once_with({'GT': '0/1'})

    @patch(
        "aliquotmaf.subcommands.vcf_to_aliquot.extractors.genotypes.GenotypeAndDepthsExtractor._extract_sanger_pindel",
        autospec=True,
        return_value=([5, 6], 11),
    )
    def test_extract_sanger_pindel(self, extract_fn):
        """
        Test dispatch Sanger Pindel
        """
        # Prepare
        caller = variant_callers.SANGER_PINDEL
        # Test
        newgt, depths = GenotypeAndDepthsExtractor._dispatch_extractor(
            1, {'GT': '0/1'}, [], caller
        )
        # Validate results
        extract_fn.assert_called_once_with({'GT': '0/1'})

    @patch(
        "aliquotmaf.subcommands.vcf_to_aliquot.extractors.genotypes.GenotypeAndDepthsExtractor._extract_svaba_somatic",
        autospec=True,
        return_value=([5, 6], 11),
    )
    def test_extract_svaba_somatic(self, extract_fn):
        """
        Test dispatch SvABA Somatic
        """
        # Prepare
        caller = variant_callers.SVABA_SOMATIC
        # Test
        newgt, depths = GenotypeAndDepthsExtractor._dispatch_extractor(
            1, {'GT': '0/1'}, [], caller
        )
        # Validate results
        extract_fn.assert_called_once_with({'GT': '0/1'})

    @patch(
        "aliquotmaf.subcommands.vcf_to_aliquot.extractors.genotypes.GenotypeAndDepthsExtractor._extract_strelka_somatic",
        autospec=True,
        return_value=([5, 6], 11),
    )
    def test_extract_strelka_somatic(self, extract_fn):
        """
        Test dispatch Strelka Somatic
        """
        # Prepare
        caller = variant_callers.STRELKA_SOMATIC
        # Test
        newgt, depths = GenotypeAndDepthsExtractor._dispatch_extractor(
            1, {'GT': '0/1'}, [], caller
        )
        # Validate results
        extract_fn.assert_called_once_with({'GT': '0/1'})

    def test__extract_mutect2_dp_set(self):
        """
        MuTect2 depth extraction algorithm
        Case where 'DP' is set properly
        """
        # Prepare
        genotype = {
            'GT': '0/1',
            'AD': (122, 45),
            'AF': 0.233,
            'DP': 167,
            'F1R2': (37, 13),
            'F2R1': (38, 10),
            'FAD': (78, 23),
            'SB': (55, 67, 22, 23),
        }
        # Test
        allele_depths, dp = GenotypeAndDepthsExtractor._extract_mutect2(genotype)
        # Expected results
        expected_allele_depths = [122, 45]
        expected_dp = 167
        # Validate results
        assert allele_depths == expected_allele_depths
        assert dp == expected_dp

    def test__extract_mutect2_dp_unset(self):
        """
        MuTect2 depth extraction algorithm
        Case where 'DP' is not set
        """
        # Prepare
        genotype = {
            'GT': '0/1',
            'AD': (122, 45),
            'AF': 0.233,
            'F1R2': (37, 13),
            'F2R1': (38, 10),
            'FAD': (78, 23),
            'SB': (55, 67, 22, 23),
        }
        # Test
        allele_depths, dp = GenotypeAndDepthsExtractor._extract_mutect2(genotype)
        # Expected results
        expected_allele_depths = [122, 45]
        expected_dp = 167
        # Validate results
        assert allele_depths == expected_allele_depths
        assert dp == expected_dp

    def test__extract_mutect2_single_allele(self):
        """
        MuTect2 depth extraction algorithm
        Case where 'DP' single allele is encountered
        """
        # Prepare
        genotype = {
            'GT': '0/0',
            'AD': 122,
            'AF': 0.0,
            'DP': 122,
            'F1R2': (37, 0),
            'F2R1': (38, 0),
            'FAD': (78, 0),
            'SB': (55, 67, 0, 0),
        }
        # Test
        allele_depths, dp = GenotypeAndDepthsExtractor._extract_mutect2(genotype)
        # Expected results
        expected_allele_depths = [122]
        expected_dp = 122
        # Validate results
        assert allele_depths == expected_allele_depths
        assert dp == expected_dp

    def test__extract_muse_single_allele(self):
        """
        MuSE depth extraction algorithm
        """
        # Prepare
        genotype = {'GT': '0/0', 'AD': 15, 'BQ': 30, 'DP': 15, 'SS': '.'}
        # Test
        allele_depths, dp = GenotypeAndDepthsExtractor._extract_muse(genotype)
        # Expected results
        expected_allele_depths = [15]
        expected_dp = 15
        # Validate results
        assert allele_depths == expected_allele_depths
        assert dp == expected_dp

    def test__extract_muse_multi_allele(self):
        """
        MuSE depth extraction algorithm
        """
        # Prepare
        genotype = {
            'GT': '0/1',
            'AD': (7, 1),
            'BQ': (33, 40),
            'DP': 8,
            'SS': '2',
        }
        # Test
        allele_depths, dp = GenotypeAndDepthsExtractor._extract_muse(genotype)
        # Expected results
        expected_allele_depths = [7, 1]
        expected_dp = 8
        # Validate results
        assert allele_depths == expected_allele_depths
        assert dp == expected_dp

    def test__extract_caveman(self):
        """
        CaVEMan depth extraction algorithm
        """
        # Prepare
        idx = 1
        genotype = {
            'GT': '0|1',
            'FAZ': 0,
            'FCZ': 0,
            'FGZ': 2,
            'FTZ': 12,
            'RAZ': 0,
            'RCZ': 0,
            'RGZ': 2,
            'RTZ': 11,
            'PM': 0.15,
        }
        alleles = ('T', 'G')
        # Test
        allele_depths, dp = GenotypeAndDepthsExtractor._extract_caveman(
            idx, genotype, alleles
        )
        # Expected results
        expected_allele_depths = [23, 4]
        expected_dp = 27
        # Validate results
        assert allele_depths == expected_allele_depths
        assert dp == expected_dp

    def test__extract_somaticsniper(self):
        """
        SomaticSniper depth extraction algorithm
        """
        # Prepare
        idx = 1
        genotype = {
            'GT': '0/0',
            'IGT': '0/0',
            'DP': 6,
            'DP4': (2, 4, 0, 0),
            'BCOUNT': (0, 0, 0, 6),
            'GQ': 45,
            'JGQ': 42,
            'VAQ': 0,
            'BQ': 33,
            'MQ': 55,
            'AMQ': 55,
            'SS': 0,
            'SSC': '.',
        }
        alleles = ('T', 'C')
        # Test
        allele_depths, dp = GenotypeAndDepthsExtractor._extract_somaticsniper(
            idx, genotype, alleles
        )
        # Expected results
        expected_allele_depths = [6, 0]
        expected_dp = 6
        # Validate results
        assert allele_depths == expected_allele_depths
        assert dp == expected_dp

    def test__extract_varscan2(self):
        """
        VarScan2 depth extraction algorithm
        """
        # Prepare
        idx = 1
        genotype = {
            'GT': '0/1',
            'GQ': '.',
            'DP': 10,
            'RD': 6,
            'AD': 4,
            'FREQ': '40%',
            'DP4': (2, 4, 0, 4),
        }
        alleles = ('G', 'A')
        # Test
        allele_depths, dp = GenotypeAndDepthsExtractor._extract_varscan2(
            idx, genotype, alleles
        )
        # Expected results
        expected_allele_depths = [6, 4]
        expected_dp = 10
        # Validate results
        assert allele_depths == expected_allele_depths
        assert dp == expected_dp

    def test__extract_pindel(self):
        """
        Pindel depth extraction algorithm
        """
        # Prepare
        genotype = {
            'GT': '0/1',
            'AD': (22, 2),
            'DP': 24,
            'SS': 2,
        }
        # Test
        allele_depths, dp = GenotypeAndDepthsExtractor._extract_pindel(genotype)
        # Expected results
        expected_allele_depths = [22, 2]
        expected_dp = 24
        # Validate results
        assert allele_depths == expected_allele_depths
        assert dp == expected_dp

    def test__extract_sanger_pindel(self):
        """
        Sanger Pindel depth extraction algorithm
        """
        # Prepare
        genotype = {
            'GT': '0/0',
            'PP': 0,
            'NP': 5,
            'PB': 0,
            'NB': 0,
            'PD': 37,
            'ND': 58,
            'PR': 37,
            'NR': 63,
            'PU': 0,
            'NU': 5,
            'FD': 98,
            'FC': 5,
        }
        # Test
        allele_depths, dp = GenotypeAndDepthsExtractor._extract_sanger_pindel(genotype)
        # Expected results
        expected_allele_depths = [95, 5]
        expected_dp = 100
        # Validate results
        assert allele_depths == expected_allele_depths
        assert dp == expected_dp

    def test__extract_svaba_somatic(self):
        """
        SvABA depth extraction algorithm
        """
        # Prepare
        genotype = {
            'GT': '0/1',
            'AD': (6,),
            'CR': 6,
            'DP': 67,
            'GQ': 2,
            'LO': 11.38,
            'LR': -1.982,
            'SR': 6,
        }
        # Test
        allele_depths, dp = GenotypeAndDepthsExtractor._extract_svaba_somatic(genotype)
        # Expected results
        expected_allele_depths = [61, 6]
        expected_dp = 67
        # Validate results
        assert allele_depths == expected_allele_depths
        assert dp == expected_dp
