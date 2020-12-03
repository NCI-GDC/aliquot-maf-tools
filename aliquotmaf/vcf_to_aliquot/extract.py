def extract(
    tumor_sample_id,
    normal_sample_id,
    tumor_idx,
    normal_idx,
    ann_cols,
    vep_key,
    record,
    is_tumor_only,
):
    """
    Extract the VCF information needed to transform into MAF.
    """
    dic = {
        "var_allele_idx": None,
        "tumor_gt": None,
        "tumor_depths": None,
        "normal_gt": None,
        "normal_depths": None,
        "location_data": None,
        "effects": None,
        "selected_effect": None,
        "variant_class": None,
    }

    # Genotypes
    var_allele_idx = Extractors.VariantAlleleIndexExtractor.extract(
        tumor_genotype=record.samples[tumor_sample_id]
    )
    tumor_gt, tumor_depths = Extractors.GenotypeAndDepthsExtractor.extract(
        var_allele_idx=var_allele_idx,
        genotype=record.samples[tumor_sample_id],
        alleles=record.alleles,
    )

    if not is_tumor_only:
        normal_gt, normal_depths = Extractors.GenotypeAndDepthsExtractor.extract(
            var_allele_idx=var_allele_idx,
            genotype=record.samples[normal_sample_id],
            alleles=record.alleles,
        )
    else:
        normal_gt, normal_depths = None, None

    # Locations
    location_data = Extractors.LocationDataExtractor.extract(
        ref_allele=record.ref,
        var_allele=record.alleles[var_allele_idx],
        position=record.pos,
        alleles=record.alleles,
    )

    # Handle effects
    effects = Extractors.EffectsExtractor.extract(
        effect_priority=self.effect_priority,
        biotype_priority=self.biotype_priority,
        effect_keys=ann_cols,
        effect_list=[urllib.parse.unquote(i).split("|") for i in record.info[vep_key]],
        var_idx=var_allele_idx,
    )

    effects, selected_effect = Extractors.SelectOneEffectExtractor.extract(
        all_effects=effects,
        effect_priority=self.effect_priority,
        biotype_priority=self.biotype_priority,
        custom_enst=self.custom_enst,
    )

    selected_effect = Extractors.PopulationFrequencyExtractor.extract(
        effect=selected_effect, var_allele=location_data["var_allele"]
    )

    # Handle variant class
    variant_class = Extractors.VariantClassExtractor.extract(
        cons=selected_effect["One_Consequence"],
        var_type=location_data["var_type"],
        inframe=location_data["inframe"],
    )

    # Make return dictionary
    dic["var_allele_idx"] = var_allele_idx
    dic["tumor_gt"] = tumor_gt
    dic["tumor_depths"] = tumor_depths
    dic["normal_gt"] = normal_gt
    dic["normal_depths"] = normal_depths
    dic["location_data"] = location_data
    dic["effects"] = format_all_effects(effects)
    dic["selected_effect"] = selected_effect
    dic["variant_class"] = variant_class
    dic["vcf_columns"] = format_vcf_columns(
        vcf_record=record, vep_key=vep_key, tumor_idx=tumor_idx, normal_idx=normal_idx,
    )
    return dic
