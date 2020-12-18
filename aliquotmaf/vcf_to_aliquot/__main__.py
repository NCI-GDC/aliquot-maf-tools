#!/usr/bin/env python3
"""
Subcommand for converting a VEP annotated VCF to a raw
aliquot MAF.
"""
import argparse
import logging
import sys
from collections import simplenamespace

from aliquotmaf.subcommands.vcf_to_aliquot.runners import GDC_1_0_0_Aliquot

logger = logging.getLogger(__name__)


def setup_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-vcf", required=True, help="Path to input VCF file")
    parser.add_argument("--output-maf", required=True, help="Path to output MAF file")

    vcf = parser.add_argument_group(title="VCF options")
    vcf.add_argument(
        "--tumor-only", action="store_true", help="Is this a tumor-only VCF?"
    )
    vcf.add_argument(
        "--tumor-vcf-id", default="TUMOR", help="Name of the tumor sample in the VCF",
    )
    vcf.add_argument(
        "--normal-vcf-id",
        default="NORMAL",
        help="Name of the normal sample in the VCF",
    )
    # TODO: Make choice from enum
    vcf.add_argument(
        "--caller-id",
        required=True,
        help="Name of the caller used to detect mutations",
    )
    vcf.add_argument(
        "--src-vcf-uuid", required=True, help="The UUID of the src VCF file"
    )

    sample = parser.add_argument_group(title="Sample Metadata")
    sample.add_argument("--case-uuid", required=True, help="Sample case UUID")
    sample.add_argument(
        "--tumor-submitter-id", required=True, help="Tumor sample aliquot submitter ID",
    )
    sample.add_argument(
        "--tumor-aliquot-uuid", required=True, help="Tumor sample aliquot UUID"
    )
    sample.add_argument("--tumor-bam-uuid", required=True, help="Tumor sample bam UUID")

    sample.add_argument(
        "--normal-submitter-id", help="Normal sample aliquot submitter ID"
    )
    sample.add_argument("--normal-aliquot-uuid", help="Normal sample aliquot UUID")
    sample.add_argument("--normal-bam-uuid", help="Normal sample bam UUID")

    sample.add_argument("--sequencer", action="append", help="The sequencer used")
    sample.add_argument(
        "--maf-center", action="append", required=True, help="The sequencing center"
    )

    anno = parser.add_argument_group(title="Annotation Resources")
    anno.add_argument(
        "--biotype-priority-file", required=True, help="Biotype priority JSON"
    )
    anno.add_argument(
        "--effect-priority-file", required=True, help="Effect priority JSON"
    )
    anno.add_argument(
        "--custom-enst", default=None, help="Optional custom ENST overrides"
    )
    anno.add_argument(
        "--dbsnp-priority-db", default=None, help="DBSNP priority sqlite database"
    )
    anno.add_argument("--reference-fasta", required=True, help="Reference fasta file")
    anno.add_argument(
        "--reference-fasta-index", required=True, help="Reference fasta fai file"
    )
    anno.add_argument(
        "--reference-context-size",
        type=int,
        default=5,
        help="Number of BP to add both upstream and "
        + "downstream from variant for reference context",
    )
    anno.add_argument(
        "--cosmic-vcf", default=None, help="Optional COSMIC VCF for annotating"
    )
    anno.add_argument(
        "--non-tcga-exac-vcf",
        default=None,
        help="Optional non-TCGA ExAC VCF for annotating and filtering",
    )
    anno.add_argument("--hotspot-tsv", default=None, help="Optional hotspot TSV")

    filt = parser.add_argument_group(title="Filtering Options")
    filt.add_argument(
        "--exac-freq-cutoff",
        default=0.001,
        type=float,
        help="Flag variants where the allele frequency in any ExAC population is great than this value as common_in_exac [0.001]",
    )
    filt.add_argument(
        "--gdc-blacklist",
        type=str,
        default=None,
        help="The file containing the blacklist tags and tumor aliquot uuids to apply them to.",
    )
    filt.add_argument(
        "--min-n-depth",
        default=7,
        type=int,
        help="Flag variants where normal depth is <= INT as ndp [7].",
    )
    filt.add_argument(
        "--gdc-pon-vcf",
        type=str,
        default=None,
        help="The tabix-indexed panel of normals VCF for applying the gdc pon filter",
    )
    filt.add_argument(
        "--nonexonic-intervals",
        type=str,
        default=None,
        help="Flag variants outside of this tabix-indexed bed file as NonExonic",
    )
    filt.add_argument(
        "--target-intervals",
        action="append",
        help="Flag variants outside of these tabix-indexed bed files as off_target. Use one or more times.",
    )
    return parser


def process_argv(argv=None) -> simplenamespace:
    pass


def run():
    # Setup annotator classes
    """
    # Initalize aliquot class
        # MafHeader: reference_fasta_index, tumor_only, tumor_aliquot_uuid, normal_aliquot_uuid
        # MafSorter: reference_fasta_index
        # MafWriter: output_maf
    for annotator in aliquot.annotators:
        annotator_instance = init_annotator(annotator)
    with Vcf() as vcf_object:
        for record in vcf_object:
            do_stuff
    # Initalize GDC_1_0_0_Aliquot
    with GDC_1_0_0_Aliquot() as aliquot, Vcf() as vcf:
        for record in vcf:
            extracted_data = record.extract(tumor_vcf_id, normal_vcf_id)
            transformed_data = record.transform(extracted_data)
            blop
    """
    pass


def main(argv=None) -> int:

    exit_code = 0
    argv = argv or sys.argv
    args = process_argv(argv)

    setup_logger(args)

    try:
        exit_code = run(args)
    except Exception as e:
        logger.exception(e)
        exit_code = 1

    # call runner
    # cls = options.func(options)
    # cls.do_work()

    return exit_code


if __name__ == '__main__':
    exit_code = 0
    try:
        exit_code = main()
    except Exception as e:
        logger.exception(e)
        exit_code = 1
    sys.exit(exit_code)

# __END__