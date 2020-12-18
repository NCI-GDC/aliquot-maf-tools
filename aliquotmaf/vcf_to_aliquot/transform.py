#!/usr/bin/env python3

from operator import itemgetter

from maflib.record import MafRecord
from maflib.validation import ValidationStringency

from aliquotmaf.vcf_to_aliquot.converters.builder import get_builder
from aliquotmaf.vcf_to_aliquot.converters.collection import InputCollection
from aliquotmaf.vcf_to_aliquot.converters.formatters import (
    format_all_effects,
    format_alleles,
    format_depths,
    format_vcf_columns,
)
from aliquotmaf.vcf_to_aliquot.filter import filter
from aliquotmaf.vcf_to_aliquot.vcf import VcfRecord

# TODO: Build input collection in own module


def transform(vcf_record: VcfRecord, data: dict, aliquot: Aliquot, line_number: int):

    # is_tumor_only, line_number=None, colset, maf_center, tumor_submitter_id, sequencer, tumor_aliquot_uuid, normal_aliquot_uuid, scheme, annotators):
    """
    Transform into maf record.
    """

    # Generic data
    collection = InputCollection()
    keys = itemgetter("selected_effect", itemgetter("Hugo_Symbol"))
    collection.add(
        column="Hugo_Symbol",
        value=data["selected_effect"].get("Hugo_Symbol"),
        default="Unknown",
    )
    collection.add(
        column="Entrez_Gene_Id", value=data["selected_effect"]["Entrez_Gene_Id"]
    )
    collection.add(column="Center", value=aliquot.maf_center)
    collection.add(column="NCBI_Build", value="GRCh38")
    collection.add(column="Chromosome", value=vcf_record.chrom)
    collection.add(column="Start_Position", value=data["location_data"]["start"])
    collection.add(column="End_Position", value=data["location_data"]["stop"])
    collection.add(column="Strand", value="+")
    collection.add(column="Variant_Classification", value=data["variant_class"])
    collection.add(column="Variant_Type", value=data["location_data"]["var_type"])
    collection.add(column="Reference_Allele", value=data["location_data"]["ref_allele"])

    for k, v in zip(
        ["Tumor_Seq_Allele1", "Tumor_Seq_Allele2"],
        format_alleles(
            genotype=data["tumor_gt"],
            alleles=data["location_data"]["alleles"],
            defaults=[
                data["location_data"]["ref_allele"],
                data["location_data"]["var_allele"],
            ],
        ),
    ):
        collection.add(column=k, value=v)

    if not aliquot.is_tumor_only:
        for k, v in zip(
            ["Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2"],
            format_alleles(
                genotype=data["normal_gt"],
                alleles=data["location_data"]["alleles"],
                defaults=[
                    data["location_data"]["ref_allele"],
                    data["location_data"]["ref_allele"],
                ],
            ),
        ):
            collection.add(column=k, value=v)
    else:
        for k in ["Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2"]:
            collection.add(column=k, value="")

    collection.add(column="dbSNP_RS", value=data["selected_effect"]["dbSNP_RS"])

    collection.add(column="Tumor_Sample_Barcode", value=aliquot.tumor_submitter_id)
    collection.add(
        column="Matched_Norm_Sample_Barcode",
        value=aliquot.normal_submitter_id,
        default="",
    )
    collection.add(column="Sequencer", value=aliquot.sequencer, default="")
    collection.add(column="Tumor_Sample_UUID", value=aliquot.tumor_aliquot_uuid)
    collection.add(
        column="Matched_Norm_Sample_UUID",
        value=aliquot.normal_aliquot_uuid,
        default="",
    )
    collection.add(column="all_effects", value=";".join(data["effects"]))

    for k, v in zip(
        ["t_depth", "t_ref_count", "t_alt_count"],
        format_depths(
            genotype=data["tumor_gt"],
            depths=data["tumor_depths"],
            var_allele_idx=data["var_allele_idx"],
            default_total_dp=0,
        ),
    ):
        collection.add(column=k, value=v)

    if not aliquot.is_tumor_only:
        for k, v in zip(
            ["n_depth", "n_ref_count", "n_alt_count"],
            format_depths(
                genotype=data["normal_gt"],
                depths=data["normal_depths"],
                var_allele_idx=data["var_allele_idx"],
            ),
        ):
            collection.add(column=k, value=v)
    else:
        for k in ["n_depth", "n_ref_count", "n_alt_count"]:
            collection.add(column=k, value=None)

    for k in data["selected_effect"]:
        if k in aliquot.colset and k not in collection._colset:
            collection.add(column=k, value=data["selected_effect"][k])

    # Set other uuids
    collection.add(column="src_vcf_id", value=aliquot.src_vcf_uuid)
    collection.add(column="tumor_bam_uuid", value=aliquot.tumor_bam_uuid)
    collection.add(column="normal_bam_uuid", value=aliquot.normal_bam_uuid)
    collection.add(column="case_id", value=aliquot.case_uuid)

    # VCF columns
    collection.add(column="FILTER", value=";".join(sorted(list(vcf_record.filter))))
    collection.add(column="vcf_region", value=data["vcf_columns"]["vcf_region"])
    collection.add(column="vcf_info", value=data["vcf_columns"]["vcf_info"])
    collection.add(column="vcf_format", value=data["vcf_columns"]["vcf_format"])
    collection.add(column="vcf_tumor_gt", value=data["vcf_columns"]["vcf_tumor_gt"])
    collection.add(
        column="vcf_normal_gt", value=data["vcf_columns"].get("vcf_normal_gt")
    )

    # Set the other columns to none
    collection.add(column="Score", value="")
    collection.add(column="BAM_File", value="")
    collection.add(column="Sequencing_Phase", value="")

    anno_set = ("dbSNP_Val_Status", "COSMIC", "CONTEXT", "Mutation_Status")
    for i in aliquot.colset - set(collection.columns()):
        if i not in anno_set:
            collection.add(column=i, value=None)
    collection.transform(aliquot.scheme)

    maf_record = MafRecord(
        line_number=line_number, stringency=ValidationStringency.Strict
    )
    for i in collection:
        maf_record += i.transformed

    # TODO: Look into using MafRecord.add(MafColumnRecord) for adding columns
    # DbSnp Validation
    dbsnp_val_status_record = aliquot.annotators["DbSnpValidation"].annotate(maf_record)
    maf_record["dbSNP_Val_Status"] = dbsnp_val_status_record

    # Cosmic
    if aliquot.annotators["cosmicId"]:
        cosmic_return = aliquot.annotators["cosmicId"].annotate(maf_record, vcf_record)
        maf_record["COSMIC"] = cosmic_return.cosmic
        if cosmic_return.dbsnp is not None:
            maf_record["dbSNP_RS"] = cosmic_return.dbsnp
    else:
        maf_record["COSMIC"] = get_builder("COSMIC", aliquot.scheme, value=None)

    # Non-TCGA Exac
    if aliquot.annotators["NonTcgaExac"]:
        non_tcga_exac_records = aliquot.annotators["NonTcgaExac"].annotate(
            maf_record, vcf_record, var_allele_idx=data["var_allele_idx"]
        )
        for pop, val in non_tcga_exac_records.items():
            maf_record[pop] = val

    # Hotspots
    if aliquot.annotators["hotspots"]:
        hotspot_record = aliquot.annotators["hotspots"].annotate(maf_record)
    else:
        hotspot_record = get_builder("hotspot", aliquot.scheme, value=None)
    maf_record['hotspot'] = hotspot_record

    # Reference context
    context_record = aliquot.annotators["reference_context"].annotate(
        maf_record, vcf_record
    )
    maf_record["CONTEXT"] = context_record

    # Mutation Status
    mutation_status_record = aliquot.annotators["mutation_status"].annotate(
        maf_record, vcf_record, aliquot.tumor_sample_id
    )
    maf_record['mutation_status'] = mutation_status_record

    # Filters
    gdc_filter_record = filter(maf_record, aliquot)
    maf_record["GDC_FILTER"] = gdc_filter_record

    return maf_record


# __END__
