"""
Tests the extractors in aliquotmaf.subcommands.vcf_to_protected.extractors.effects
"""

import pytest

from aliquotmaf.subcommands.vcf_to_aliquot.extractors.effects import (
    EffectsExtractor,
    EffectsExtractor_102,
    SelectOneEffectExtractor,
)

# VariantAlleleIndexExtractor -> GenotypeAndDepthsExtractor -> LocationDataExtractor -> EffectsExtractor -> SelectOneEffectExtractor -> PopulationFrequencyExtractor -> VariantClassExtractor

# VEP Columns in current impl
ANNO_COLUMNS = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|RefSeq|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|ENTREZ|EVIDENCE".split(
    "|"
)

ANNO_COLUMNS_102 = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|RefSeq|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS".split(
    "|"
)


def create_basic_effects(allele, allele_num):
    lst = []

    effects = ["A", "B", "C"]
    biotypes = ["A", "B", "C"]
    for i in effects:
        for j in biotypes:
            rec = []
            for col in ANNO_COLUMNS:
                if col == "Consequence":
                    rec.append(i)
                elif col == "BIOTYPE":
                    rec.append(j)
                elif col == "Allele":
                    rec.append(allele)
                elif col == "ALLELE_NUM":
                    rec.append(allele_num)
                else:
                    rec.append("")

            lst.append(rec)
    return lst


@pytest.mark.parametrize(
    "effect_priority, biotype_priority, var_idx, enst, expected",
    [
        (
            {"A": 1, "B": 2, "C": 3, "": 10},
            {"A": 1, "B": 2, "C": 3},
            1,
            None,
            ("A", "A"),
        ),
        (
            {"A": 1, "B": 2, "C": 3, "": 10},
            {"A": 2, "B": 1, "C": 3},
            1,
            None,
            ("A", "B"),
        ),
        (
            {"A": 1, "B": 2, "C": 3, "": 10},
            {"A": 2, "B": 3, "C": 1},
            1,
            None,
            ("A", "C"),
        ),
        (
            {"A": 3, "B": 1, "C": 3, "": 10},
            {"A": 1, "B": 2, "C": 3},
            1,
            None,
            ("B", "A"),
        ),
        (
            {"A": 1, "B": 1, "C": 3, "": 10},
            {"A": 2, "B": 1, "C": 3},
            1,
            None,
            ("A", "B"),
        ),
    ],
)
def test_select_one_effect_extractor(
    effect_priority, biotype_priority, var_idx, enst, expected
):
    """
    Tests the extraction of single effect based on the priority data. This is
    only an extremely simple test, and more work needs to be done to test all
    the various nuances of the implementation.
    """
    elist = create_basic_effects("A", var_idx)
    res = EffectsExtractor.extract(
        effect_priority, biotype_priority, ANNO_COLUMNS, elist, var_idx
    )
    all_effects, selected = SelectOneEffectExtractor.extract(
        res, effect_priority, biotype_priority, custom_enst=enst
    )
    assert (selected["Consequence"], selected["BIOTYPE"]) == expected


@pytest.mark.parametrize(
    "effect_priority, biotype_priority, var_idx, enst, expected",
    [
        (
            {"A": 1, "B": 2, "C": 3, "": 10},
            {"A": 1, "B": 2, "C": 3},
            1,
            None,
            ("A", "A"),
        ),
        (
            {"A": 1, "B": 2, "C": 3, "": 10},
            {"A": 2, "B": 1, "C": 3},
            1,
            None,
            ("A", "B"),
        ),
        (
            {"A": 1, "B": 2, "C": 3, "": 10},
            {"A": 2, "B": 3, "C": 1},
            1,
            None,
            ("A", "C"),
        ),
        (
            {"A": 3, "B": 1, "C": 3, "": 10},
            {"A": 1, "B": 2, "C": 3},
            1,
            None,
            ("B", "A"),
        ),
        (
            {"A": 1, "B": 1, "C": 3, "": 10},
            {"A": 2, "B": 1, "C": 3},
            1,
            None,
            ("A", "B"),
        ),
    ],
)
def test_select_one_effect_extractor_vep_102(
    effect_priority, biotype_priority, var_idx, enst, expected
):
    """
    Tests the extraction of single effect based on the priority data. This is
    only an extremely simple test, and more work needs to be done to test all
    the various nuances of the implementation.
    """
    elist = create_basic_effects("A", var_idx)
    res = EffectsExtractor_102.extract(
        effect_priority, ANNO_COLUMNS_102, elist, var_idx
    )
    all_effects, selected = SelectOneEffectExtractor.extract(
        res, effect_priority, biotype_priority, custom_enst=enst
    )
    assert (selected["Consequence"], selected["BIOTYPE"]) == expected
