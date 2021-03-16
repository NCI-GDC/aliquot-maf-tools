"""
Extractor class for genomic region data.
"""
from typing import List, NamedTuple

from aliquotmaf.vcf_to_aliquot.extractors.base import Extractor


class LocationNT(NamedTuple):
    alleles: List[str]
    var_allele: str
    ref_allele: str
    start: int
    stop: int
    var_type: str
    inframe: bool


class LocationDataExtractor(Extractor):
    """
    Extracts the start, stop, variant type, and inframe values from the provided
    data. Mutates alleles based on normalizing the ref and alt alleles.
    """

    NP_TYPE = {1: "SNP", 2: "DNP", 3: "TNP"}

    @classmethod
    def extract(cls, ref_allele, var_allele, position, alleles) -> LocationNT:
        ref_length, var_length = len(ref_allele), len(var_allele)

        # Remove any prefixed reference bps from all alleles, using "=" for simple indels
        while (
            ref_allele
            and var_allele
            and ref_allele[0] == var_allele[0]
            and ref_allele != var_allele
        ):

            ref_allele = ref_allele[1:] if len(ref_allele) > 1 else "-"
            var_allele = var_allele[1:] if len(var_allele) > 1 else "-"
            alleles = [i[1:] if len(i) > 1 else "-" for i in alleles]
            ref_length -= 1
            var_length -= 1
            position += 1

        # position variables
        start, stop, var_type, inframe = None, None, None, None

        # Handle SNPs, DNPs, TNPs, or anything larger
        if ref_length == var_length:
            start, stop = position, position + var_length - 1
            var_type = cls.NP_TYPE.get(var_length, "ONP")

        # Handle all indels, including those complex ones which contain substitutions
        elif ref_length != var_length:
            # Insertions
            if ref_length < var_length:
                start = position - 1 if ref_allele == "-" else position
                stop = position if ref_allele == "-" else position + ref_length - 1
                var_type = "INS"
            # Deletions
            else:
                start, stop = position, position + ref_length - 1
                var_type = "DEL"

            inframe = abs(ref_length - var_length) % 3 == 0

        rdic = {
            "alleles": alleles,
            "var_allele": var_allele,
            "ref_allele": ref_allele,
            "start": start,
            "stop": stop,
            "var_type": var_type,
            "inframe": inframe,
        }
        return LocationNT(**rdic)
