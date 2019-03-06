"""
Tests the extractors in aliquotmaf.subcommands.vcf_to_aliquot.extractors.population_frequency
"""
import pytest

from aliquotmaf.subcommands.vcf_to_aliquot.extractors.population_frequency import PopulationFrequencyExtractor

# VariantAlleleIndexExtractor -> GenotypeAndDepthsExtractor -> LocationDataExtractor -> EffectsExtractor -> SelectOneEffectExtractor -> PopulationFrequencyExtractor -> VariantClassExtractor

SUBPOPS = ['GMAF', 'AFR_MAF', 'AMR_MAF', 'ASN_MAF', 'EAS_MAF',
  'EUR_MAF', 'SAS_MAF', 'AA_MAF', 'EA_MAF', 'ExAC_AF_Adj',
  'ExAC_AF', 'ExAC_AF_AFR', 'ExAC_AF_AMR', 'ExAC_AF_EAS',
  'ExAC_AF_FIN', 'ExAC_AF_NFE', 'ExAC_AF_OTH', 'ExAC_AF_SAS']

@pytest.mark.parametrize("effect, var_allele, expected", [
    ({'GMAF': 'A:0.1', 'ExAC_AF_NFE': 'A:0.2'}, 'A', {'GMAF': '0.1', 'ExAC_AF_NFE': '0.2'}),
    ({'GMAF': 'A:0.1;T:0.0003', 'ExAC_AF_NFE': 'A:0.2'}, 'A', {'GMAF': '0.1', 'ExAC_AF_NFE': '0.2'}),
    ({'GMAF': 'T:0.0003', 'ExAC_AF_NFE': 'A:0.2'}, 'A', {'ExAC_AF_NFE': '0.2'}),
    ({'GMAF': 'A:0.1;T:0.0003', 'ExAC_AF_NFE': 'A:0.2'}, 'T', {'GMAF': '0.0003'}),
])
def test_population_frequency_extractor(effect, var_allele, expected):
    """
    Tests the extraction of population frequencies. 
    """
    raw = PopulationFrequencyExtractor.extract(effect, var_allele) 
    res = {k : raw[k] for k in raw if raw[k] is not None} 
    assert res == expected
