"""
Applies the "common in non-cancer" gnomAD filter.
"""

from __future__ import absolute_import

from .filter_base import Filter


class FilterGnomAD(Filter):
    def __init__(self, cutoff):
        super().__init__(name="CommonInGnomAD")
        self.tags = ["common_in_gnomAD"]
        self.cutoff = cutoff
        self.maxAF_field = "gnomAD_non_cancer_MAX_AF_adj"
        self.logger.info(
            "Using gnomAD non-cancer MAX_AF_adj frequency cutoff of {0}".format(cutoff)
        )

    @classmethod
    def setup(cls, cutoff):
        curr = cls(cutoff)
        return curr

    def filter(self, maf_record):
        freq = maf_record[self.maxAF_field].value
        if freq is not None and freq > self.cutoff:
            return True
        return False

    def shutdown(self) -> None:
        pass
