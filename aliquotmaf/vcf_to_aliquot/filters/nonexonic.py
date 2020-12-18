"""
Applies the NonExonic filter for regions outside the gencode intervals provided.
"""
from pysam import TabixFile, asBed

from aliquotmaf.filters.filter_base import Filter


class NonExonic(Filter):
    def __init__(self, source):
        super().__init__(name="NonExonic", source=source)
        self.tags = ["NonExonic"]
        self.f = None
        self.logger.info("Using genode exon interval file {0}".format(source))

    @classmethod
    def setup(cls, args):
        source = args.nonexonic_intervals
        if not source:
            return None
        curr = cls(source)
        curr.f = TabixFile(curr.source, parser=asBed())
        return curr

    def filter(self, maf_record):
        flag = True
        vcf_region = maf_record["vcf_region"].value.split(":")
        pos = int(vcf_region[1])
        region = "{0}:{1}-{2}".format(vcf_region[0], pos, maf_record["End_Position"])
        try:
            for record in self.f.fetch(region=region):
                flag = False
                break
        except ValueError:
            pass
        return flag

    def shutdown(self):
        self.f.close()
