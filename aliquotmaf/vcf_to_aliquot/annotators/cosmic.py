"""
Annotates the COSMIC ID and mutates the dbSNP_RS if necessary.
"""

from typing import NamedTuple, Optional

import pysam
from maflib.column import MafColumnRecord

from aliquotmaf.vcf_to_aliquot.annotators.annotator import Annotator
from aliquotmaf.vcf_to_aliquot.converters.builder import get_builder


class CosmicId(Annotator):
    def __init__(self, scheme, source):
        super().__init__(name="CosmicID", source=source, scheme=scheme)
        self.f = None

    @classmethod
    def setup(cls, scheme, args):
        curr = cls(scheme, args.cosmic_vcf)
        curr.f = pysam.VariantFile(curr.source)
        return curr

    def annotate(self, maf_record, vcf_record, var_allele_idx=1) -> MafColumnRecord:
        region = "{0}:{1}-{2}".format(
            vcf_record.chrom, vcf_record.pos, vcf_record.pos + 1
        )
        alt = vcf_record.alleles[var_allele_idx]
        cosmic_ids = []
        for record in self.f.fetch(region=region):
            try:
                if (
                    vcf_record.pos == record.pos
                    and vcf_record.ref == record.ref
                    and alt == record.alts[0]
                ):
                    cosmic_ids.append(record.id)
            except TypeError:
                # Weirdly formatted COSMIC variants
                pass

        # NOTE: Overwrite `dbSNP_rs` if cosmic IDs found and dbSNP reports `novel`
        if cosmic_ids:
            cosmic_record = get_builder(
                "COSMIC", self.scheme, value=";".join(sorted(list(set(cosmic_ids))))
            )
        else:
            cosmic_record = get_builder("COSMIC", self.scheme, value=None)

        return cosmic_record

    def shutdown(self):
        self.f.close()
