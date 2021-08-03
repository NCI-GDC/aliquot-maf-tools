"""
Implements the gnomAD annotation using the gnomAD VCF.
"""

from __future__ import absolute_import

import os
from collections import OrderedDict
from itertools import repeat
from sys import path

import pysam

from aliquotmaf.converters.builder import get_builder

from .annotator import Annotator

GNOMAD_SRC_TO_MAF = OrderedDict(
    AF_non_cancer_eas="gnomAD_non_cancer_EAS_AF",
    AF_non_cancer_afr="gnomAD_non_cancer_AFR_AF",
    AF_non_cancer_ami="gnomAD_non_cancer_AMI_AF",
    AF_non_cancer_mid="gnomAD_non_cancer_MID_AF",
    AF_non_cancer_sas="gnomAD_non_cancer_SAS_AF",
    AF_non_cancer_nfe="gnomAD_non_cancer_NFE_AF",
    AF_non_cancer="gnomAD_non_cancer_AF",
    AF_non_cancer_amr="gnomAD_non_cancer_AMR_AF",
    AF_non_cancer_oth="gnomAD_non_cancer_OTH_AF",
    AF_non_cancer_asj="gnomAD_non_cancer_ASJ_AF",
    AF_non_cancer_fin="gnomAD_non_cancer_FIN_AF",
    MAX_AF_non_cancer_adj="gnomAD_non_cancer_MAX_AF_adj",
    POP_MAX_non_cancer_adj="gnomAD_non_cancer_MAX_AF_POPS_adj",
)
GNOMAD_SOURCE_COLUMNS = GNOMAD_SRC_TO_MAF.keys()
GNOMAD_MAF_COLUMNS = GNOMAD_SRC_TO_MAF.values()


class GnomAD_VCF(Annotator):
    def __init__(self, scheme, source):
        super().__init__(name="GnomAD", scheme=scheme, source=source)
        self.f = None

    @classmethod
    def setup(cls, scheme, source):
        curr = cls(scheme, source)
        curr.f = pysam.VariantFile(curr.source)
        return curr

    def annotate(self, maf_record, vcf_record, var_allele_idx=1):
        """
        Annotate each variant with AF records from GnomAD
        """
        region = "{0}:{1}-{2}".format(vcf_record.chrom, vcf_record.pos, vcf_record.stop)
        alt = vcf_record.alleles[var_allele_idx]
        found = False

        for record in self.f.fetch(region=region):
            if (
                vcf_record.pos == record.pos
                and vcf_record.ref == record.ref
                and alt in record.alts
            ):
                for source_col, maf_col in GNOMAD_SRC_TO_MAF.items():
                    found = True
                    value = record.info.get(source_col)
                    default = ''

                    if source_col == "POP_MAX_non_cancer_adj" and value is not None:
                        value = list(value)
                        default = []

                    elif type(value) == tuple:
                        value = value[0]

                    maf_record[maf_col] = get_builder(
                        maf_col, self.scheme, value=value, default=default
                    )

                break

        if found:
            return maf_record
        else:
            for source_col, maf_col in GNOMAD_SRC_TO_MAF.items():
                default = ''

                if source_col == "POP_MAX_non_cancer_adj":
                    default = []

                maf_record[maf_col] = get_builder(maf_col, self.scheme, value=default)
            return maf_record

    def shutdown(self):
        """Close the pysam VCF object"""
        self.f.close()
