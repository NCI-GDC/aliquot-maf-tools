"""
Applies the normal depth filter
"""
from aliquotmaf.vcf_to_aliquot.filters.filter_base import Filter


class NormalDepth(Filter):
    def __init__(self, cutoff, is_tumor_only: bool = False):
        super().__init__(name="NormalDepth")
        self.tags = ["ndp"] if not is_tumor_only else []
        self.cutoff = cutoff

    @classmethod
    def setup(cls, args):
        curr = cls(args.min_n_depth, args.is_tumor_only)
        return curr

    def filter(self, maf_record):
        ndp = maf_record["n_depth"].value
        return ndp is not None and ndp <= self.cutoff

    def shutdown(self):
        pass
