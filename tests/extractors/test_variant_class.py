"""
Tests the extractors in aliquotmaf.subcommands.vcf_to_aliquot.extractors.variant_class
"""
import pytest

from aliquotmaf.subcommands.vcf_to_aliquot.extractors.variant_class import (
    VariantClassExtractor,
)


@pytest.mark.parametrize(
    "cons, var_type, inframe, expected",
    [
        ("splice_donor_variant", "", "", "Splice_Site"),
        ("exon_loss_variant", "", "", "Splice_Site"),
        ("stop_gained", "", "", "Nonsense_Mutation"),
        ("frameshift_variant", "DEL", "", "Frame_Shift_Del"),
        ("protein_altering_variant", "DEL", False, "Frame_Shift_Del"),
        ("frameshift_variant", "INS", "", "Frame_Shift_Ins"),
        ("protein_altering_variant", "INS", False, "Frame_Shift_Ins"),
        ("stop_lost", "", "", "Nonstop_Mutation"),
        ("start_lost", "", "", "Translation_Start_Site"),
        ("inframe_insertion", "INS", "", "In_Frame_Ins"),
        ("protein_altering_variant", "INS", True, "In_Frame_Ins"),
        ("inframe_deletion", "DEL", "", "In_Frame_Del"),
        ("protein_altering_variant", "DEL", True, "In_Frame_Del"),
        ("missense_variant", "", "", "Missense_Mutation"),
        ("intron_variant", "", "", "Intron"),
        ("splice_region_variant", "", "", "Splice_Region"),
        ("synonymous_variant", "", "", "Silent"),
        ("non_coding_exon_variant", "", "", "RNA"),
        ("5_prime_UTR_variant", "", "", "5'UTR"),
        ("3_prime_UTR_variant", "", "", "3'UTR"),
        ("intergenic_variant", "", "", "IGR"),
        ("upstream_gene_variant", "", "", "5'Flank"),
        ("downstream_gene_variant", "", "", "3'Flank"),
        ("random", "", "", "Targeted_Region"),
    ],
)
def test_variant_class_extractor(cons, var_type, inframe, expected):
    """
    Tests the extraction of MAF variant class. 
    """
    res = VariantClassExtractor.extract(cons, var_type, inframe)
    assert res == expected
