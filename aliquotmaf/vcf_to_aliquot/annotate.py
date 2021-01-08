#!/usr/bin/env python3

from maflib.record import MafRecord

from aliquotmaf.vcf_to_aliquot.aliquot.base import Aliquot
from aliquotmaf.vcf_to_aliquot.converters.builder import get_builder


def annotate(
    aliquot: Aliquot, data: dict, maf_record: MafRecord, vcf_record
) -> MafRecord:
    # TODO: Return list of new mafrecords to add
    # TODO: Business logic for records updating current records
    # TODO: Look into using MafRecord.add(MafColumnRecord) for adding columns
    """
    annotation_cols = {
        name: annotator.annotate(data, maf_record, vcf_record)
        for name, annotator in aliquot.annotators.items()

    }
    """

    # DbSnp Validation
    dbsnp_val_status_record = aliquot.annotators["DbSnpValidation"].annotate(maf_record)
    maf_record["dbSNP_Val_Status"] = dbsnp_val_status_record

    # Cosmic
    # Set 'novel' dbSNP annotation to None if cosmic returns results
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

    return maf_record


# __END__
