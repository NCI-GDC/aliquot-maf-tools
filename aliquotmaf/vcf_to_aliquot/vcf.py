#!/usr/bin/env python3

from types import SimpleNamespace
from typing import Generator, List, Optional

import pysam
from maflib.record import MafRecord

from aliquotmaf.vcf_to_aliquot import extract, transform

DI = SimpleNamespace(pysam=pysam)


class VcfRecord:
    def __init__(
        self, record, index: int = None,
    ):
        self.record = record
        self.index = index


class VcfFile:
    def __init__(
        self,
        variant_file: str,
        tumor_vcf_id: str = None,
        normal_vcf_id: str = None,
        is_tumor_only: bool = False,
        di=DI,
    ):
        self.variant_file = variant_file
        self.tumor_vcf_id = tumor_vcf_id
        self.normal_vcf_id = normal_vcf_id
        self.is_tumor_only = is_tumor_only

        self.pysam_file = None

        self.tumor_idx = self.assert_sample_id_in_header(tumor_vcf_id)
        self.normal_idx = self.assert_sample_id_in_header(normal_vcf_id)
        self.annotation_columns = self.extract_annotation_from_header()

        self.vep_key = "CSQ"
        self.di = di

    def __enter__(self):
        self.pysam_file = self.di.pysam.VariantFile(self.variant_file)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.pysam_file.close()

    def __iter__(self) -> Generator[VcfRecord, None, None]:
        counter = 0
        for record in self.pysam_file.fetch():
            counter += 1
            yield VcfRecord(
                record=record, index=counter,
            )

    def assert_sample_id_in_header(self, vcf_id) -> Optional[int]:
        """
        Asserts that a given sample is in the VCF header and returns the index.

        :param vcf_object: ``~pysam.VariantFile`` instance
        :param sample: the sample's name
        :return: sample column index
        """
        idx = None

        try:
            slist = list(self.pysam_file.header.samples)
            idx = slist.index(vcf_id) + 9  # FIXME: Why plus 9 here, header?
        except ValueError:
            pass
        return idx

    def extract_annotation_from_header(self, vep_key="CSQ") -> List[str]:
        """
        Extract the VEP annotation columns from the VCF header.

        :param vcf_object: ``~pysam.VariantFile`` instance
        :return: annotation column list
        :return: the key in the INFO data from the VCF
        """
        # annotation cols
        ann_cols_format = []

        # Loop over header records
        for record in self.pysam_file.header.records:
            if record.type == "INFO":
                iname = record.get("ID")
                if iname and str(iname) == vep_key:  # This could be a problem
                    vep_key = iname  # following records tested against previous
                    anno_line = re.search(r'Format: (\S+)"$', record["Description"])
                    raw_ann_cols_format = anno_line.group(1).split("|")
                    for ann in raw_ann_cols_format:
                        if ann.startswith("ExAC"):
                            ann_cols_format.append(self.fix_exac(ann))
                        else:
                            ann_cols_format.append(ann)

        return ann_cols_format

    def fix_exac(self, ann):
        """
        Convert the to the "OLD" version of ExAC headers from the old plugin.
        before VEP included it built in.

        :param ann: the annotation key
        :return: the converted annotation key
        """
        lookup = {
            "ExAC_MAF": "ExAC_AF",
            "ExAC_AFR_MAF": "ExAC_AF_AFR",
            "ExAC_AMR_MAF": "ExAC_AF_AMR",
            "ExAC_EAS_MAF": "ExAC_AF_EAS",
            "ExAC_FIN_MAF": "ExAC_AF_FIN",
            "ExAC_NFE_MAF": "ExAC_AF_NFE",
            "ExAC_OTH_MAF": "ExAC_AF_OTH",
            "ExAC_SAS_MAF": "ExAC_AF_SAS",
            "ExAC_Adj_MAF": "ExAC_AF_Adj",
        }
        return lookup[ann]


# __END__
