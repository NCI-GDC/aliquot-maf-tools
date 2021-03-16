"""
Tests for the ``aliquotmaf.subcommands.vcf_to_aliquot`` subcommands.
"""
import uuid

import pytest

from aliquotmaf.vcf_to_aliquot.__main__ import main


def test_validation_tumor_only():
    fake_uuid = str(uuid.uuid4())
    resource_commands = [
        "--biotype-priority-file",
        "./fake.json",
        "--effect-priority-file",
        "./fake.json",
        "--reference-fasta",
        "./fake.fa",
        "--reference-fasta-index",
        "./fake.fa.fai",
    ]

    main_commands = [
        "VcfToAliquotMaf",
        "--input-vcf",
        "./fake.vcf.gz",
        "--output-maf",
        "./fake.maf.gz",
        "gdc-1.0.0-aliquot",
        "--caller-id",
        "mutect",
        "--src-vcf-uuid",
        fake_uuid,
        "--case-uuid",
        fake_uuid,
        "--tumor-submitter-id",
        "FAKE-TUMOR",
        "--tumor-aliquot-uuid",
        fake_uuid,
        "--tumor-bam-uuid",
        fake_uuid,
        "--maf-center",
        "FAKE",
    ] + resource_commands

    with pytest.raises(ValueError):
        main(main_commands)

    with pytest.raises(FileNotFoundError):
        main(main_commands + ["--tumor-only"])
