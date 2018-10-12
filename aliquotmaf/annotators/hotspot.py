"""
Implements the hotspots annotation. 
"""
from __future__ import absolute_import

from .annotator import Annotator
from aliquotmaf.converters.builder import get_builder

class Hotspot(Annotator):
    def __init__(self, source, scheme, data):
        super().__init__(name='Hotspot', source=source, scheme=scheme)
        self.data = data 

    @classmethod
    def setup(cls, scheme, source):
        # load the hotspots

        hsdic = {}
        head = []
        count = 0
        with open(source, 'rt') as fh:
            for line in fh:
                if not head:
                    head = line.rstrip('\r\n').lower().split('\t')
                    assert all([i in head for i in ['hugo_symbol', 'change', 'type']]), \
                        self.logger.error("Unexpected header {0} found!".format(head))
                else:
                    dat = dict(zip(head, line.rstrip('\r\n').split('\t')))
                    if dat['hugo_symbol'] not in hsdic: 
                        hsdic[dat['hugo_symbol']] = {}
                    hsdic[dat['hugo_symbol']][dat['change']] = dat['type']
                    count += 1 
        curr = cls(source, scheme, hsdic)
        curr.logger.info("Loaded {0} hotspots".format(count))
        return curr

    def annotate(self, maf_record):
        gene = maf_record['Hugo_Symbol'].value
        mval = "N" 
        if gene in self.data:
            hgvsp = None if not maf_record['HGVSp_Short'].value else \
                maf_record['HGVSp_Short'].value.lstrip('p.')
            if hgvsp and hgvsp in self.data[gene]:
                self.logger.info(hgvsp)
                mval = "Y"
        maf_record['hotspot'] = get_builder("hotspot", self.scheme, value=mval)
        return maf_record

    def shutdown(self):
        pass 
