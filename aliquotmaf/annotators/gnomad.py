"""
Implements the gnomAD annotation.

*** DEPRECATED AND UNUSED ***
This implementation had performance and scaling issues and was
abandoned in favor of the vcf implementation
"""

from __future__ import absolute_import

import os
from collections import OrderedDict
from itertools import repeat

import pandas as pd

from aliquotmaf.converters.builder import get_builder

from .annotator import Annotator

CHROM_LIST = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
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


def empty_gnomad_record():
    return pd.Series(dict(zip(GNOMAD_SOURCE_COLUMNS, repeat(""))))


class GnomAD(Annotator):
    def __init__(self, scheme, ref_prefix):
        super().__init__(name="GnomAD", scheme=scheme)
        self.file_template = ref_prefix + "{}.feather"
        self.chrom = None
        self.df = None

    @classmethod
    def setup(cls, scheme, ref_prefix):
        """
        Prepare annotation sources

        ref_prefix: prefix used to locate reference files it will be prepended
            to strings in the form of 'chr1.feather' in order to construct
            the path to the reference files
        """
        curr = cls(scheme, ref_prefix)
        missing = [
            chr
            for chr in CHROM_LIST
            if not os.path.exists(curr.file_template.format(chr))
        ]
        if missing:
            curr.logger.warning(
                "Missing reference files for chromosomes {}".format(missing)
            )
        else:
            curr.logger.info(
                "Found GnomAD reference files using prefix {}".format(ref_prefix)
            )
        return curr

    def get_gnomad_record(self, chrom, pos, ref, alt):
        """
        Manage in-memory data and retrieve annotations
        """

        if chrom is not self.chrom:
            del self.df
            self.df = self.load_chrom(chrom)
            self.chrom = chrom
        return self.lookup_variant(pos, ref, alt)

    def lookup_variant(self, pos, ref, alt) -> pd.Series:
        """
        lookup variant record in currently loaded dataframe
        """
        if pos not in self.df.index:
            # position has no gnomad annotation -> return empty record
            return empty_gnomad_record()
        vdf = self.df.loc[pos]
        if isinstance(vdf, pd.Series):
            # single record at that position
            target = f"{ref}|{alt}"
            if vdf["ref_alt"] == target:
                # record matches variant
                return vdf[list(GNOMAD_SOURCE_COLUMNS)].fillna("")
            else:
                # does not match variant
                return empty_gnomad_record()
        elif isinstance(vdf, pd.DataFrame):
            # multiple records at that position
            vdf = vdf.reset_index()
            target = f"{ref}|{alt}"
            # do as little work as possible to find a match
            for idx, value in vdf["ref_alt"].items():
                if value == target:
                    return vdf.loc[idx][list(GNOMAD_SOURCE_COLUMNS)].fillna("")
            # failed to find match
            return empty_gnomad_record()

    def load_chrom(self, chrom):
        """
        Load chromosome into memory
        """
        if chrom in CHROM_LIST:
            filepath = self.file_template.format(chrom)
            return pd.read_feather(filepath).set_index("pos")
        else:
            raise KeyError("Unrecognized contig encountered: {}".format(chrom))

    def annotate(self, maf_record, vcf_record):
        """
        Annotate each variant with AF records from GnomAD
        """
        grec = self.get_gnomad_record(
            vcf_record.chrom, vcf_record.pos, vcf_record.ref, vcf_record.alts[0]
        )

        for source_col, maf_col in GNOMAD_SRC_TO_MAF.items():
            value = str(getattr(grec, source_col))
            default = ""
            if source_col == "POP_MAX_non_cancer_adj":
                value = value.split(",")
                default = []
            maf_record[maf_col] = get_builder(
                maf_col, self.scheme, value=value, default=default
            )
        return maf_record

    def shutdown(self):
        """
        This annotator has no open files to close.
        """
        pass
