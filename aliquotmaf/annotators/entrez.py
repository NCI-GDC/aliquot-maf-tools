"""
Implements the Entrez annotation.
"""

from __future__ import absolute_import

from json import load

from maflib.schemes import MafScheme

from aliquotmaf.converters.builder import get_builder

from .annotator import Annotator

MAF_SYMBOL: str = "Hugo_Symbol"
MAF_FEATURE: str = "Feature"


class Entrez(Annotator):
    def __init__(self, scheme, source):
        super().__init__(name="Entrez", source=source, scheme=scheme)
        self.gencode: dict
        self.ncbi: dict

    @classmethod
    def setup(cls, scheme: MafScheme, source: str) -> "Entrez":
        """
        Load annotation data
        """
        curr = cls(scheme, source)
        with open(curr.source, "r") as fh:
            tmp = load(fh)
        curr.gencode = tmp["GENCODE"]
        curr.ncbi = tmp["NCBI"]
        curr.logger.info(
            "Loaded {} GENCODE to ENTREZ mappings".format(len(curr.gencode))
        )
        curr.logger.info("Loaded {} NCBI to ENTREZ mappings".format(len(curr.gencode)))
        return curr

    def annotate(self, maf_record):
        """
        Annotate provided record with Entrez gene ID if possible
        """
        # field names like "SYMBOL" and "Feature" should come from schema
        # default value in case no match is found
        entrez_id = [0]
        # try to find match based on HGNC Symbol first
        if MAF_SYMBOL in maf_record:
            ncbi = maf_record[MAF_SYMBOL].value
            entrez_id = self.ncbi.get(ncbi, [0])
        # otherwise try to find match based on Ensembl Transcript ID
        elif MAF_FEATURE in maf_record:
            gencode = maf_record[MAF_FEATURE].value
            entrez_id = self.gencode.get(gencode, [0])
        # Replicating old VEP entrez filter behavior, if multiple entrez IDs
        # are available just select the first one in the list
        maf_record["Entrez_Gene_Id"] = get_builder(
            "Entrez_Gene_Id", self.scheme, value=str(entrez_id[0]), default="0"
        )
        return maf_record

    def shutdown(self):
        """
        Annotator end of life actions
        """
        # nothing to do here
        pass
