#!/usr/bin/env python3

import urllib
from typing import List, NamedTuple, Optional

from aliquotmaf.utils.utils import load_enst, load_json
from aliquotmaf.vcf_to_aliquot import extractors
from aliquotmaf.vcf_to_aliquot.aliquot.base import Aliquot
from aliquotmaf.vcf_to_aliquot.converters.formatters import (
    format_all_effects,
    format_vcf_columns,
)
from aliquotmaf.vcf_to_aliquot.vcf import VcfRecord


class ExtractedDataNT(NamedTuple):
    var_allele_idx: int
    tumor_gt: str
    tumor_depths: List[int]
    normal_gt: str
    normal_depths: List[int]
    location_data: dict
    effects: List[str]
    selected_effect: str
    variant_class: str


def extract(
    vcf_record: VcfRecord,
    aliquot: Aliquot,
    tumor_idx: int,
    biotype_priority_file: str,
    effect_priority_file: str,
    annotation_columns: List[str],
    normal_idx: Optional[int] = None,
    custom_enst: Optional[str] = None,
    vep_key: str = "CSQ",
    _extractors=extractors,
) -> ExtractedDataNT:
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
    record = vcf_record.record
    tumor_sample_id = aliquot.tumor_sample_id
    normal_sample_id = aliquot.normal_sample_id
    is_tumor_only = aliquot.is_tumor_only

    # Genotypes
    # TODO: Combine
    var_allele_idx = _extractors.VariantAlleleIndexExtractor.extract(
        tumor_genotype=record.samples[tumor_sample_id]
    )
    # TODO: Improve returns
    tumor_gt, tumor_depths = _extractors.GenotypeAndDepthsExtractor.extract(
        var_allele_idx=var_allele_idx,
        genotype=record.samples[tumor_sample_id],
        alleles=record.alleles,
    )

    if not is_tumor_only:
        # TODO: Improve returns
        normal_gt, normal_depths = _extractors.GenotypeAndDepthsExtractor.extract(
            var_allele_idx=var_allele_idx,
            genotype=record.samples[normal_sample_id],
            alleles=record.alleles,
        )
    else:
        normal_gt, normal_depths = None, None

    # Locations
    location_data = _extractors.LocationDataExtractor.extract(
        ref_allele=record.ref,
        var_allele=record.alleles[var_allele_idx],
        position=record.pos,
        alleles=record.alleles,
    )
    location_data = location_data._as_dict()

    effect_priority = load_json(effect_priority_file)
    biotype_priority = load_json(biotype_priority_file)
    # Handle effects
    effects = _extractors.EffectsExtractor.extract(
        effect_priority=effect_priority,
        biotype_priority=biotype_priority,
        effect_keys=annotation_columns,
        effect_list=[urllib.parse.unquote(i).split("|") for i in record.info[vep_key]],
        var_idx=var_allele_idx,
    )

    # Combine next two
    # TODO: Improve returns
    effects, selected_effect = _extractors.SelectOneEffectExtractor.extract(
        all_effects=effects,
        effect_priority=effect_priority,
        biotype_priority=biotype_priority,
        custom_enst=load_enst(custom_enst),
    )

    selected_effect = _extractors.PopulationFrequencyExtractor.extract(
        effect=selected_effect, var_allele=location_data["var_allele"]
    )

    # Handle variant class
    variant_class = _extractors.VariantClassExtractor.extract(
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
    dic["effects"] = format_all_effects(effects)  # FIXME: Call elsewhere
    dic["selected_effect"] = selected_effect
    dic["variant_class"] = variant_class
    dic["vcf_columns"] = format_vcf_columns(
        vcf_record=record, vep_key=vep_key, tumor_idx=tumor_idx, normal_idx=normal_idx,
    )  # FIXME: Call elsewhere
    return dic
