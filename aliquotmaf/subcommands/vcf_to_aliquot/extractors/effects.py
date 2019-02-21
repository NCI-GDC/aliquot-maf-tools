"""Module for extractor classes that help extract and format the effect data from
some input source, mostly engineered for VEP annotated VCF files.

* EffectsExtractor         Extract the VEP effects into a list of dictionaries
* SelectOneEffectExtractor Select a single transcript effect based on priorities
"""
import re

from aliquotmaf.subcommands.vcf_to_aliquot.extractors import Extractor

class EffectsExtractor(Extractor):
    """A `~maf_converter_lib.extractor.Extractor` class that takes the VEP
    annotated data in formats into a list of dictionaries.
    """
    # Long to short amino-acids translation
    AA3TO1= {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Asx': 'B',
             'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Glx': 'Z', 'Gly': 'G',
             'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M',
             'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
             'Tyr': 'Y', 'Val': 'V', 'Xxx': 'X', 'Ter': '*'}

    @classmethod
    def extract(cls, effect_priority, biotype_priority,
              effect_keys, effect_list, var_idx):
        """
        Entry point for parsing VEP effects.

        :param effect_priority: A `dict` of effect types and their priority
        :param biotype_priority: A `dict` of biotypes and their priority
        :param effect_keys: `list` of the keys corresponding to the values in
                            effect_list
        :param effect_list: `list` of lists of effect data
        :param var_idx: the variant allele index
        :returns: a `list` of effect `dict`
        """
        # The list of parsed effects to return
        all_effects = []

        # Loop over each effect
        for edat in effect_list:
            effect = {}
            for i,v in enumerate(effect_keys):
                try:
                    if edat[i]: effect[v] = edat[i].replace('&', ';')
                    else: effect[v] = None
                except: effect[v] = None

            # Skip effects on other ALT alleles. If ALLELE_NUM is undefined (e.g.
            # for INFO:SVTYPE), don't skip any
            if effect['ALLELE_NUM'] and int(effect['ALLELE_NUM']) != var_idx:
                continue

            # Remove transcript ID from HGVS codon/protein changes,
            # to make it easier on the eye
            if effect['HGVSc']: effect['HGVSc'] = re.sub(r'^.*:', '', effect["HGVSc"])
            if effect['HGVSp']: effect['HGVSp'] = re.sub(r'^.*:', '', effect["HGVSp"])


            # Remove the predefined HGVSc code in HGVSp, if found
            if effect['HGVSp'] and effect['HGVSp'].startswith('c.'):
                effect['HGVSp'] = re.sub(r'^.*\((p\.\S+)\)', r'\1', effect['HGVSp'])

            # If there are several consequences listed for a transcript,
            # choose the most severe one
            effect['One_Consequence'] = sorted(effect['Consequence'].split(';'),
                                               key = lambda x: effect_priority.get(x, 20))[0] \
                                               if effect['Consequence'] \
                                               else 'intergenic_variant'

            # Create a shorter HGVS protein format using 1-letter codes
            if effect['HGVSp']:
                hgvs_p_short = effect['HGVSp']
                for aa in cls.AA3TO1: hgvs_p_short = hgvs_p_short.replace(aa, cls.AA3TO1[aa])
                effect['HGVSp_Short'] = hgvs_p_short

            # Fix HGVSp_Short, CDS_position, and Protein_position for splice
            # acceptor/donor variants
            if re.search(r'^(splice_acceptor_variant|splice_donor_variant)$',
                        effect['One_Consequence']) and effect['HGVSc']:
                c_pos = re.search(r'^c.(\d+)', effect['HGVSc'])
                if c_pos is not None:
                    c_pos = 1 if int(c_pos.group(1)) < 1 else int(c_pos.group(1))
                    p_pos = '{0:.0f}'.format((c_pos + c_pos % 3) / 3.0)
                    effect['HGVSp_Short']      = "p.X" + p_pos + "_splice"
                    effect['CDS_position']     = re.sub(r'^-(/\d+)$', str(c_pos) + r'\1',
                                                        effect['CDS_position'])
                    effect['Protein_position'] = re.sub(r'^-(/\d+)$', str(p_pos) + r'\1',
                                                    effect['Protein_position'])

            # Fix HGVSp_Short for Silent mutations, so it mentions the amino-acid and position
            if 'HGVSp_Short' in effect and effect['HGVSp_Short'] == "p.=":
                p_pos = re.search(r'^(\d+)(-\d+)?/\d+$', effect['Protein_position']).group(1)
                aa    = effect['Amino_acids']
                effect['HGVSp_Short'] = 'p.{0}{1}{0}'.format(aa, p_pos)

            # Copy VEP data into MAF fields that don't share the same identifier
            effect['Transcript_ID'] = effect['Feature']
            effect['Exon_Number']   = effect['EXON']
            effect['Hugo_Symbol']   = effect['SYMBOL']

            # GDC Add ENTREZ ID
            if effect['ENTREZ']:
                # ::NOTE:: Mapping is ambiguous in some cases. I am just shifting off one id
                # if there are many.
                effect['Entrez_Gene_Id'] = effect['ENTREZ'].split(';')[0]
            else: effect['Entrez_Gene_Id'] = None

            # If VEP couldn't find this variant in dbSNP/etc., we'll say it's "novel"
            if effect['Existing_variation']:
                # ::NOTE:: If seen in a DB other than dbSNP, this field will remain blank
                effect['dbSNP_RS'] = ';'.join([i for i in \
                    effect['Existing_variation'].split(';') if re.search(r'^rs\d+$', i)])
                if effect['dbSNP_RS'] == '': effect['dbSNP_RS'] = None
            else:
                effect['dbSNP_RS'] = 'novel'

            # Transcript_Length isn't separately reported, but can be parsed out
            # from cDNA_position
            if effect['cDNA_position']:
                tlen = re.search(r'\/(\d+)$', effect['cDNA_position'])
                effect['Transcript_Length'] = tlen.group(1) if tlen is not None else '0'
            else:
                effect['Transcript_Length'] = '0'

            # Append the current effect to the all_effects list
            all_effects.append(effect)

        return all_effects

class SelectOneEffectExtractor(Extractor):
    """A `~maf_converter_lib.extractor.Extractor` class that takes the priority
    data and selects a single effect.
    """
    @classmethod
    def extract(cls, all_effects, effect_priority, biotype_priority,
              custom_enst=None):

        maf_effect = None

        # Sort effects first by transcript biotype, then by severity,
        # and then by longest transcript
        all_effects = sorted(
            all_effects,
            key = lambda x: (biotype_priority.get(x['BIOTYPE'], 20),
                             effect_priority.get(x['One_Consequence'], 20),
                             -1 * int(x['Transcript_Length'])))

        # Find the highest priority effect with a gene symbol (usually the first one)
        try: maf_effect_with_gene_name = list(filter(lambda x: x['SYMBOL'], all_effects))[0]
        except IndexError: maf_effect_with_gene_name = None
        maf_gene = maf_effect_with_gene_name['SYMBOL'] \
            if maf_effect_with_gene_name is not None else ''

        # If the gene has user-defined custom isoform overrides, choose that instead
        if custom_enst:
            try: maf_effect = list(filter(lambda x: x['SYMBOL'] and x['SYMBOL'] == maf_gene \
                                          and x['Transcript_ID'] and x['Transcript_ID'] \
                                          in custom_enst, all_effects))[0]
            except IndexError: maf_effect = None

        if not maf_effect:
            # Find the effect on the canonical transcript of that highest priority gene
            try: maf_effect = list(filter(lambda x: x['SYMBOL'] and x['SYMBOL'] == maf_gene \
                                and x['CANONICAL'] and x['CANONICAL'] == 'YES',
                                all_effects))[0]
            except IndexError: maf_effect = None

            # If that gene has no canonical transcript tagged, choose the highest priority
            # canonical effect on any gene
            if not maf_effect:
                try: maf_effect = list(filter(lambda x: x['CANONICAL'] == 'YES', all_effects))[0]
                except IndexError: maf_effect = None

            # If none of the effects are tagged as canonical, then just report top priority effect
            if not maf_effect:
                maf_effect = all_effects[0]

        return all_effects, maf_effect
