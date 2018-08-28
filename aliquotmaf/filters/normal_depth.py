"""
Applies the normal depth filter 
"""
from __future__ import absolute_import

from .filter_base import Filter

class NormalDepth(Filter):
    def __init__(self, cutoff):
        super().__init__(name='NormalDepth')
        self.tags = ['ndp']
        self.cutoff = cutoff
        self.logger.info("Using normal depth cutoff of {0}".format(cutoff))

    @classmethod
    def setup(cls, cutoff):
        curr = cls(cutoff)
        return curr

    def filter(self, maf_record):
        ndp = maf_record['n_depth'].value
        return ndp is not None and ndp <= self.cutoff

    def shutdown(self):
        pass
