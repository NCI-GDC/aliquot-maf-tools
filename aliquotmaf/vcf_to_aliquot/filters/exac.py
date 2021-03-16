"""
Applies the common in non-tcga ExAC filter.
"""

from aliquotmaf.vcf_to_aliquot.filters.filter_base import Filter


class ExAC(Filter):
    def __init__(self, cutoff):
        super().__init__(name="CommonInExAC")
        self.tags = ["common_in_exac"]
        self.cutoff = cutoff
        self.subpops = [
            "nontcga_ExAC_AF_Adj",
            "nontcga_ExAC_AF",
            "nontcga_ExAC_AF_AFR",
            "nontcga_ExAC_AF_AMR",
            "nontcga_ExAC_AF_EAS",
            "nontcga_ExAC_AF_FIN",
            "nontcga_ExAC_AF_NFE",
            "nontcga_ExAC_AF_OTH",
            "nontcga_ExAC_AF_SAS",
        ]

    @classmethod
    def setup(cls, args):
        curr = cls(args.exac_freq_cutoff)
        return curr

    def filter(self, maf_record):
        flag = False
        for subpop in self.subpops:
            freq = maf_record[subpop].value
            if freq is not None and freq > self.cutoff:
                flag = True
                break
        return flag

    def shutdown(self):
        pass
