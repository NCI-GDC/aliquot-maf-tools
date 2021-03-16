#!/usr/bin/env python3

from typing import Callable

from maflib.record import MafRecord
from maflib.validation import ValidationStringency

from aliquotmaf.vcf_to_aliquot.aliquot import Aliquot
from aliquotmaf.vcf_to_aliquot.annotate import annotate
from aliquotmaf.vcf_to_aliquot.converters.collection import InputCollection
from aliquotmaf.vcf_to_aliquot.converters.formatters import (
    format_alleles,
    format_depths,
)
from aliquotmaf.vcf_to_aliquot.filter_maf import filter_maf
from aliquotmaf.vcf_to_aliquot.vcf import VcfRecord


def initalize_collection(vcf_record, data, aliquot):
    collection = InputCollection()
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
    collection.add(
        column="FILTER", value=";".join(sorted(list(vcf_record.record.filter)))
    )
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

    return collection


def transform(
    vcf_record: VcfRecord,
    data: dict,
    aliquot: Aliquot,
    line_number: int,
    _annotate: Callable = annotate,
    _filter_maf: Callable = filter_maf,
    _initalize_collection: Callable = initalize_collection,
) -> MafRecord:
    """
    Converts VCF Record to MAF Record.
    """

    record = vcf_record.record
    # Generic data
    collection = _initalize_collection(vcf_record, data, aliquot)
    collection.transform(aliquot.scheme)

    maf_record = MafRecord(
        line_number=line_number, stringency=ValidationStringency.Strict
    )
    for i in collection:
        maf_record += i.transformed

    # Annotate
    # TODO: Following
    """
    annotation_columns = _annotate()
    for col in annotation_columns:
        maf_record.add(col)
    """
    maf_record = _annotate(aliquot, data, maf_record, vcf_record=record)

    # Filters
    gdc_filter_record = _filter_maf(maf_record, aliquot)
    maf_record["GDC_FILTER"] = gdc_filter_record

    return maf_record


# __END__
