"""
Applies the GDC PON filter. We don't care about alleles, just positions. 
"""
from __future__ import absolute_import

from .filter_base import Filter

from pysam import VariantFile

class GdcPon(Filter):
    def __init__(self, source):
        super().__init__(name='GDCPON', source=source)
        self.tags = ['gdc_pon']
        self.f = None 
        self.logger.info("Using panel of normal VCF {0}".format(source))

    @classmethod
    def setup(cls, source):
        curr = cls(source)
        curr.f = VariantFile(curr.source)
        return curr

    def filter(self, maf_record):
        vcf_region = maf_record['vcf_region'].value.split(':')
        pos = int(vcf_region[1])
        region = '{0}:{1}-{2}'.format(vcf_region[0], pos, pos + 1)
        for record in self.f.fetch(region=region):
            if record.pos == pos:
                return True
        return False

    def shutdown(self):
        self.f.close()
