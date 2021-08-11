"""
Tests for the ``aliquotmaf.subcommands.vcf_to_aliquot`` subcommands.
"""
import uuid

import pytest

from aliquotmaf.__main__ import main


def test_validation_tumor_only():
    fake_uuid = str(uuid.uuid4())
    resource_commands = [
        "--biotype_priority_file",
        "./fake.json",
        "--effect_priority_file",
        "./fake.json",
        "--reference_fasta",
        "./fake.fa",
        "--reference_fasta_index",
        "./fake.fa.fai",
    ]

    main_commands = [
        "VcfToAliquotMaf",
        "--input_vcf",
        "./fake.vcf.gz",
        "--output_maf",
        "./fake.maf.gz",
        "gdc-2.0.0-aliquot",
        "--caller_id",
        "mutect",
        "--src_vcf_uuid",
        fake_uuid,
        "--case_uuid",
        fake_uuid,
        "--tumor_submitter_id",
        "FAKE-TUMOR",
        "--tumor_aliquot_uuid",
        fake_uuid,
        "--tumor_bam_uuid",
        fake_uuid,
        "--maf_center",
        "FAKE",
    ] + resource_commands

    with pytest.raises(ValueError):
        main(main_commands)

    with pytest.raises(FileNotFoundError):
        main(main_commands + ["--tumor_only"])
