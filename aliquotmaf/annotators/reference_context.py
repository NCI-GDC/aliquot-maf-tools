"""
Implements the reference context annotation.
"""
from __future__ import absolute_import

import pysam

from .annotator import Annotator
from aliquotmaf.converters.builder import get_builder

class ReferenceContext(Annotator):
    def __init__(self, source, scheme, context_size=5):
        super().__init__(name='ReferenceContext', source=source, scheme=scheme)
        self.fa = None
        self.context_size = context_size

    @classmethod
    def setup(cls, scheme, source, context_size=5):
        curr = cls(source, scheme, context_size)
        curr.fa = pysam.FastaFile(curr.source)
        return curr

    def annotate(self, maf_record, vcf_record, strip_chr=False):
        # Add reference context
        if strip_chr:
            region = '{0}:{1}-{2}'.format(
                vcf_record.chrom.replace('chr', '') if vcf_record.chrom != 'chrM' else 'MT',
                max(1, vcf_record.pos - self.context_size),
                vcf_record.stop + self.context_size
            )
        else:
            region = '{0}:{1}-{2}'.format(
                vcf_record.chrom,
                max(1, vcf_record.pos - self.context_size),
                vcf_record.stop + self.context_size
            )
        maf_record['CONTEXT'] = get_builder("CONTEXT", self.scheme, value=self.fa.fetch(region=region))
        return maf_record


    def shutdown(self):
        self.fa.close()
