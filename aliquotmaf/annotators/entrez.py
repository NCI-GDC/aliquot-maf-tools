"""
Implements the Entrez annotation.
"""

from __future__ import absolute_import

from collections import namedtuple
from json import load
from typing import TextIO

from aliquotmaf.converters.builder import get_builder

from .annotator import Annotator

# maybe expand this to a full schema at some point?
VCF_COLS = namedtuple('VCF_COLS', ['SYMBOL', 'Feature'])
Vcf = VCF_COLS('SYMBOL', 'Feature')


class Entrez(Annotator):
    def __init__(self, scheme, source):
        super().__init__(name="Entrez", source=source, scheme=scheme)

    @classmethod
    def setup(cls, scheme, source):
        """
        Load annotation data
        """
        curr = cls(scheme, source)
        with open(curr.source, 'r') as fh:
            tmp = load(fh)
        curr.gencode = tmp['GENCODE']
        curr.ncbi = tmp['NCBI']
        return curr

    def annotate(self, maf_record, vcf_record):
        """
        Annotate provided record with Entrez gene ID if possible
        """
        # field names like "SYMBOL" and "Feature" should come from schema
        # default value in case no match is found
        entrez_id = 0
        # try to find match based on HGNC Symbol first
        if Vcf.SYMBOL in vcf_record.info:
            ncbi = vcf_record.info[Vcf.SYMBOL]
            entrez_id = self.ncbi.get(ncbi, 0)
        # otherwise try to find match based on Ensembl Transcript ID
        elif Vcf.Feature in vcf_record.info:
            gencode = vcf_record.info[Vcf.Feature]
            entrez_id = self.gencode.get(gencode, 0)
        # Replicating old VEP entrez filter behavior, if multiple entrez IDs
        # are available just select the first one in the list
        maf_record["Entrez_Gene_Id"] = get_builder(
            "Entrez_Gene_Id", self.scheme, value=str(entrez_id[0]), default=0
        )
        return maf_record

    def shutdown(self):
        """
        Annotator end of life actions
        """
        # nothing to do here
        pass
