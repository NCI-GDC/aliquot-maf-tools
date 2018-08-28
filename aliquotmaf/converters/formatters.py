"""
A module containing formatting utility functions needed when converting from a 
VCF to a MAF.
"""

def format_alleles(genotype, alleles, defaults=None):
    """
    Formats the alleles from a genotype record.
    """
    al1, al2 = '', ''
    if defaults:
        al1 = defaults[0]
        al2 = defaults[1]

    if 'GT' in genotype and genotype['GT']:
        als = genotype['GT']
        if len(als) == 1:
            als = sorted(list(als) + [0])
        al1 = alleles[als[0]]
        al2 = alleles[als[1]]
    return [al1, al2]

def format_depths(genotype, depths, var_allele_idx, default_total_dp=None):
    """
    Format the variant allele depths based on the selected variant allele.
    """
    _dp = default_total_dp
    if 'DP' in genotype and \
    genotype['DP'] is not None and \
    genotype['DP'] != '.':
        _dp = genotype['DP']

    elif depths: _dp = sum(depths)

    ref_dp, alt_dp = None, None
    if depths:
        ref_dp = depths[0]
        alt_dp = depths[var_allele_idx]
    return [_dp, ref_dp, alt_dp]

def format_all_effects(all_effects):
    """
    Formats the all effects column.
    """
    # format all effects
    fmt_all_effects = []
    for eff in all_effects:
        gene_name      = eff['Hugo_Symbol'] if eff['Hugo_Symbol'] else ''
        effect_type    = eff['One_Consequence']
        protein_change = eff['HGVSp_Short'] if eff.get('HGVSp_Short') else ''
        transcript_id  = eff['Transcript_ID'] if eff.get('Transcript_ID') else ''
        # We don't want to have extra commas when there are multiple RefSeq IDs
        refseq_ids     = eff['RefSeq'].replace(';', '&') if eff.get('RefSeq') else ''
        cds_change     = eff['HGVSc'] if eff.get('HGVSc') else ''
        impact         = eff['IMPACT'] if eff.get('IMPACT') else ''
        canonical      = eff['CANONICAL'] if eff.get('CANONICAL') else ''
        sift           = eff['SIFT'].replace(';', '&') if eff.get('SIFT') else ''
        polyphen       = eff['PolyPhen'].replace(';', '&') if eff.get('PolyPhen') else ''
        strand         = eff['STRAND'] if eff.get('STRAND') else ''
        if gene_name and effect_type and transcript_id:
            curr_all_effect = ','.join(
                map(str, [gene_name, effect_type, protein_change, transcript_id, refseq_ids,
                cds_change, impact, canonical, sift, polyphen, strand]))
            fmt_all_effects.append(curr_all_effect)
    return fmt_all_effects

def format_vcf_columns(vcf_record, vep_key, tumor_idx, normal_idx=None):
    """
    Formats the VCF columns that are stored in the MAF record.
    """
    dic = {
        "vcf_region": None,
        "vcf_info": None,
        "vcf_format": None,
        "vcf_tumor_gt": None,
        "FILTER": None
    }
    # Convert pysam object to string representation
    cols = str(vcf_record).rstrip('\r\n').split('\t')
    dic['vcf_region'] = '{0}:{1}:{2}:{3}:{4}'.format(
        vcf_record.chrom,
        vcf_record.pos,
        vcf_record.id if vcf_record.id else '.',
        vcf_record.ref,
        ','.join(list(vcf_record.alts))
    )

    dic['vcf_info'] = ";".join([i for i in cols[7].split(';') \
        if not i.startswith(vep_key + '=') and \
        not i.startswith('ANN=') and not i.startswith('CSQ=')])

    dic['vcf_format'] = cols[8]
    dic['vcf_tumor_gt'] = cols[tumor_idx]
    if normal_idx is not None:
        dic['vcf_normal_gt'] = cols[normal_idx]

    dic["FILTER"] = cols[6]
    return dic
