"""
Applies the common in ExAC filter.
"""
from __future__ import absolute_import

from .filter_base import Filter

class ExAC(Filter):
    def __init__(self, cutoff):
        super().__init__(name='CommonInExAC')
        self.tags = ['common_in_exac']
        self.cutoff = cutoff
        self.subpops = ['ExAC_AF_Adj', 'ExAC_AF', 'ExAC_AF_AFR', 'ExAC_AF_AMR',
                        'ExAC_AF_EAS', 'ExAC_AF_FIN', 'ExAC_AF_NFE', 'ExAC_AF_OTH', 
                        'ExAC_AF_SAS']
        self.logger.info("Using ExAC frequency cutoff of {0}".format(cutoff))

    @classmethod
    def setup(cls, cutoff):
        curr = cls(cutoff)
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
