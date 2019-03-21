"""
Annotates the COSMIC ID and mutates the dbSNP_RS if necessary.
"""
from __future__ import absolute_import

import pysam

from .annotator import Annotator

from aliquotmaf.converters.builder import get_builder

class CosmicID(Annotator):
    def __init__(self, scheme, source):
        super().__init__(name='CosmicID', source=source, scheme=scheme)
        self.f = None

    @classmethod
    def setup(cls, scheme, source):
        curr = cls(scheme, source)
        curr.f = pysam.VariantFile(curr.source)
        return curr

    def annotate(self, maf_record, vcf_record, var_allele_idx=1):
        region = '{0}:{1}-{2}'.format(vcf_record.chrom, vcf_record.pos, vcf_record.pos + 1)
        alt = vcf_record.alleles[var_allele_idx]
        cosmic_ids = []
        for record in self.f.fetch(region=region):
            try:
                if vcf_record.pos == record.pos and \
                   vcf_record.ref == record.ref and \
                   alt == record.alts[0]:
                    cosmic_ids.append(record.id)
            except TypeError:
                # Weirdly formatted COSMIC variants
                pass

        if cosmic_ids:
            if maf_record['dbSNP_RS'].value == ['novel']:
                maf_record['dbSNP_RS'] = get_builder("dbSNP_RS", self.scheme, value=None)
            maf_record['COSMIC'] = get_builder(
                "COSMIC", self.scheme, 
                value=';'.join(sorted(list(set(cosmic_ids)))))
        else:
            maf_record['COSMIC'] = get_builder("COSMIC", self.scheme, value=None) 
        
        return maf_record 


    def shutdown(self):
        self.f.close()
