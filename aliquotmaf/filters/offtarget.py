"""
Applies the off_target filter for regions outside the provided intervals. 
"""
from __future__ import absolute_import

from .filter_base import Filter

from pysam import TabixFile, asBed


class OffTarget(Filter):
    def __init__(self, source):
        super().__init__(name="OffTarget", source=source)
        self.tags = ["off_target"]
        self.fs = []
        self.logger.info("Using interval files {0}".format(", ".join(source)))

    @classmethod
    def setup(cls, source):
        curr = cls(source)
        curr.fs = [TabixFile(i, parser=asBed()) for i in curr.source]
        return curr

    def filter(self, maf_record):
        flag = True
        vcf_region = maf_record["vcf_region"].value.split(":")
        pos = int(vcf_region[1])
        region = "{0}:{1}-{2}".format(vcf_region[0], pos, maf_record["End_Position"])
        for source in self.fs:
            try:
                for record in source.fetch(region=region):
                    return False
            except ValueError:
                pass
        return flag

    def shutdown(self):
        for i in self.fs:
            i.close()
