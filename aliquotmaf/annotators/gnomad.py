"""
Implements the gnomAD annotation.
"""

from __future__ import absolute_import

import os
from itertools import repeat
from sys import path

import pandas as pd

from aliquotmaf.converters.builder import get_builder

from .annotator import Annotator

GNOMAD_SOURCE_COLUMNS = (
    'AF_non_cancer_eas',
    'AF_non_cancer_afr',
    'AF_non_cancer_ami',
    'AF_non_cancer_mid',
    'AF_non_cancer_sas',
    'AF_non_cancer_nfe',
    'AF_non_cancer',
    'AF_non_cancer_amr',
    'AF_non_cancer_oth',
    'AF_non_cancer_asj',
    'AF_non_cancer_fin',
    'MAX_AF_non_cancer_adj',
    'POP_MAX_non_cancer_adj',
)

GNOMAD_MAF_COLUMNS = (
    "gnomAD_non_cancer_EAS_AF",
    "gnomAD_non_cancer_AFR_AF",
    "gnomAD_non_cancer_AMI_AF",
    "gnomAD_non_cancer_MID_AF",
    "gnomAD_non_cancer_SAS_AF",
    "gnomAD_non_cancer_NFE_AF",
    "gnomAD_non_cancer_AF",
    "gnomAD_non_cancer_AMR_AF",
    "gnomAD_non_cancer_OTH_AF",
    "gnomAD_non_cancer_ASJ_AF",
    "gnomAD_non_cancer_FIN_AF",
    "gnomAD_non_cancer_MAX_AF_adj",
    "gnomAD_non_cancer_MAX_AF_POPS_adj",
)

GNOMAD_SRC_TO_MAF = dict(zip(GNOMAD_SOURCE_COLUMNS, GNOMAD_MAF_COLUMNS))


class GnomAD(Annotator):
    def __init__(self, scheme, refpath, refpattern):
        super().__init__(name="GnomAD", scheme=scheme)
        self.path = refpath
        self.file_template = refpattern
        self.chrom = None
        self.chrom_list = {
            'chr1',
            'chr2',
            'chr3',
            'chr4',
            'chr5',
            'chr6',
            'chr7',
            'chr8',
            'chr9',
            'chr10',
            'chr11',
            'chr12',
            'chr13',
            'chr14',
            'chr15',
            'chr16',
            'chr17',
            'chr18',
            'chr19',
            'chr20',
            'chr21',
            'chr22',
            'chrX',
            'chrY',
            'chrM',
        }
        self.df = None

    @classmethod
    def setup(cls, scheme, refpath, refpattern):
        """
        Prepare annotation sources

        refpath should point to directory containing annotation files
        refpattern should be a template with which to create the file names
            by substituting chromosome names in for {}
            it should look like:  reference.{}.feather
        """
        curr = cls(scheme, refpath, refpattern)
        return curr

    def get_gnomad_record(self, chrom, pos, ref, alt):
        """
        Manage in-memory data and retrieve annotations
        """

        if chrom is not self.chrom:
            del self.df
            self.df = self.load_chrom(chrom)
        return self.lookup_variant(pos, ref, alt)

    def lookup_variant(self, pos, ref, alt) -> pd.Series:
        """
        lookup variant record in currently loaded dataframe
        """
        if pos not in self.df.index:
            # return empty records
            return pd.Series(dict(zip(GNOMAD_SOURCE_COLUMNS, repeat(""))))
        vdf = self.df.loc[pos]
        if type(vdf) is pd.Series:
            # single record at that position
            return vdf[list(GNOMAD_SOURCE_COLUMNS)].fillna("")
        elif type(vdf) is pd.DataFrame:
            # multiple records at that position
            vdf = vdf.reset_index()
            target = f'{ref}|{alt}'
            for idx, value in vdf['ref_alt'].iteritems():
                if value == target:
                    print(idx)
                    return vdf.loc[idx][list(GNOMAD_SOURCE_COLUMNS)].fillna("")

    def load_chrom(self, chrom):
        """
        Load chromosome into memory
        """
        if chrom in self.chrom_list:
            filepath = os.path.join(self.path, self.file_template.format(chrom))
            return pd.read_feather(filepath).set_index('pos')
        else:
            raise KeyError('Unrecognized contig encountered: {}'.format(chrom))

    def annotate(self, maf_record, vcf_record):
        """
        Annotate each variant with AF records from GnomAD
        """
        grec = self.get_gnomad_record(
            vcf_record.chrom, vcf_record.pos, vcf_record.ref, vcf_record.alts[0]
        )

        for source_col in list(GNOMAD_SOURCE_COLUMNS):
            maf_col = GNOMAD_SRC_TO_MAF[source_col]
            value = str(getattr(grec, source_col))
            default = ''
            if source_col == "POP_MAX_non_cancer_adj":
                value = value.split(',')
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
