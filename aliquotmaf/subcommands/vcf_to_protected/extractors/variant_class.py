"""
Extractor class for variant classification.
"""
import re

from aliquotmaf.subcommands.vcf_to_protected.extractors import Extractor

class VariantClassExtractor(Extractor):
    """
    Maps the selected effect to the MAF variant class options.
    """
    @classmethod
    def extract(cls, cons, var_type, inframe):

        # Splice_Site
        if re.search(r'^(splice_acceptor_variant|splice_donor_variant|transcript_ablation|exon_loss_variant)$', cons):
            return "Splice_Site"

        # stop_gained == Nonsense_Mutation
        if cons == 'stop_gained': return 'Nonsense_Mutation'

        # Frame_Shift_Del
        if (cons == 'frameshift_variant' or \
           (cons == 'protein_altering_variant' \
           and not inframe )) and var_type == 'DEL':
            return 'Frame_Shift_Del'

        # Frame_Shift_Ins
        if (cons == 'frameshift_variant' or \
           (cons == 'protein_altering_variant' and not inframe)) \
           and var_type == 'INS':
            return 'Frame_Shift_Ins'

        # stop_lost == Nonstop_Mutation
        if cons == 'stop_lost': return 'Nonstop_Mutation'

        # Translation_Start_Site
        if re.search(r'^(initiator_codon_variant|start_lost)$', cons):
            return 'Translation_Start_Site'

        # In_Frame_Ins
        if re.search(r'^(inframe_insertion|disruptive_inframe_insertion)$', cons) or \
          (cons == 'protein_altering_variant' and inframe and var_type == 'INS'):
            return 'In_Frame_Ins'

        # In_Frame_Del
        if re.search(r'^(inframe_deletion|disruptive_inframe_deletion)$', cons) or \
          (cons == 'protein_altering_variant' and inframe and var_type == 'DEL'):
            return 'In_Frame_Del'

        # Missense_Mutation
        if re.search(r'^(missense_variant|coding_sequence_variant|conservative_missense_variant|rare_amino_acid_variant)$', cons):
            return 'Missense_Mutation'

        # Introns
        if re.search(r'^(transcript_amplification|intron_variant|INTRAGENIC|intragenic_variant)$',
          cons):
            return 'Intron'

        # Slice_Region
        if cons == 'splice_region_variant':
            return 'Splice_Region'

        # Silent
        if re.search(r'^(incomplete_terminal_codon_variant|synonymous_variant|stop_retained_variant|NMD_transcript_variant)$', cons):
            return 'Silent'

        # RNA
        if re.search(r'^(mature_miRNA_variant|exon_variant|non_coding_exon_variant|non_coding_transcript_exon_variant|non_coding_transcript_variant|nc_transcript_variant)$', cons):
            return "RNA"

        # 5'UTR
        if re.search(r'^(5_prime_UTR_variant|5_prime_UTR_premature_start_codon_gain_variant)$', cons):
            return "5'UTR"

        # 3'UTR
        if cons == '3_prime_UTR_variant':
            return "3'UTR"

        # IGR
        if re.search(r'^(TF_binding_site_variant|regulatory_region_variant|regulatory_region|intergenic_variant|intergenic_region)$',
          cons):
            return "IGR"

        # 5'Flank
        if cons == 'upstream_gene_variant': return "5'Flank"

        # 3'Flank
        if cons == 'downstream_gene_variant': return "3'Flank"

        # Annotate everything else simply as a targeted region
        # TFBS_ablation, TFBS_amplification,regulatory_region_ablation,
        # regulatory_region_amplification,
        # feature_elongation, feature_truncation
        return "Targeted_Region"
