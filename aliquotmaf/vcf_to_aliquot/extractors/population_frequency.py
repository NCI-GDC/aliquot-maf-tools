"""
Extractor class for population frequency.
"""
import re
from typing import NamedTuple, Optional

from aliquotmaf.vcf_to_aliquot.extractors.base import Extractor


class PopulationFrequencyNT(NamedTuple):
    effect: Optional[str]


SUBPOPULATIONS = (
    "GMAF",
    "AFR_MAF",
    "AMR_MAF",
    "ASN_MAF",
    "EAS_MAF",
    "EUR_MAF",
    "SAS_MAF",
    "AA_MAF",
    "EA_MAF",
    "ExAC_AF_Adj",
    "ExAC_AF",
    "ExAC_AF_AFR",
    "ExAC_AF_AMR",
    "ExAC_AF_EAS",
    "ExAC_AF_FIN",
    "ExAC_AF_NFE",
    "ExAC_AF_OTH",
    "ExAC_AF_SAS",
)


class PopulationFrequencyExtractor(Extractor):
    @classmethod
    def extract(
        cls, effect, var_allele, subpops=SUBPOPULATIONS
    ) -> PopulationFrequencyNT:
        # ::NOTE:: I now correctly match alleles, in the perl version my regex didn't anchor the
        # search to the beginning of the string.

        for subpop in subpops:
            if effect.get(subpop):
                als = list(
                    filter(
                        lambda x: x.startswith(var_allele + ":"),
                        effect[subpop].split(";"),
                    )
                )
                if als:
                    val = als.pop().split(":")[1]
                    effect[subpop] = val
                else:
                    # set to None
                    effect[subpop] = None
            else:
                effect[subpop] = None
        return PopulationFrequencyNT(effect=effect)
