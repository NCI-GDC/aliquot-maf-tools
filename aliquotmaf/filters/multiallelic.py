"""
Applies the multiallelic filter.
"""

from __future__ import absolute_import

from .filter_base import Filter


class Multiallelic(Filter):
    def __init__(self):
        super().__init__(name="Multiallelic")
        self.tags = ["multiallelic"]
        self.logger.info("Loading Multialellic filter")

    @classmethod
    def setup(cls):
        curr = cls()
        return curr

    def filter(self, maf_record):
        vcf_region = maf_record["vcf_region"].value.split(":")
        alleles = set(vcf_region[3].split(",") + vcf_region[4].split(","))
        return len(alleles) > 2

    def shutdown(self):
        pass
