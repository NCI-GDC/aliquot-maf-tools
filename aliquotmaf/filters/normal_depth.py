"""
Applies the normal depth filter
"""

from __future__ import absolute_import

from .filter_base import Filter


class NormalDepth(Filter):
    def __init__(self, cutoff):
        super().__init__(name="NormalDepth")
        self.tags = ["ndp"]
        self.cutoff = cutoff
        if cutoff is not None:
            self.logger.info("Using normal depth cutoff of {0}".format(cutoff))
        else:
            self.logger.info("Normal depth cutoff not in use")

    @classmethod
    def setup(cls, cutoff):
        curr = cls(cutoff)
        return curr

    def filter(self, maf_record):
        if self.cutoff is None:
            return False
        ndp = maf_record["n_depth"].value
        return ndp is not None and ndp <= self.cutoff

    def shutdown(self):
        pass
