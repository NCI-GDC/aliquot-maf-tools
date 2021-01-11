#!/usr/bin/env python3
"""
Implements the reference context annotation.
"""
import pysam

from aliquotmaf.vcf_to_aliquot.annotators.annotator import Annotator
from aliquotmaf.vcf_to_aliquot.converters.builder import get_builder


class ReferenceContext(Annotator):
    def __init__(self, source, scheme, context_size=5):
        super().__init__(name="ReferenceContext", source=source, scheme=scheme)
        self.fa = None
        self.context_size = context_size

    @classmethod
    def setup(cls, scheme, args, default_context=5):
        curr = cls(
            args.reference_context,
            scheme,
            getattr(args, "context_size", default_context),
        )
        curr.fa = pysam.FastaFile(curr.source)
        return curr

    def annotate(self, maf_record, vcf_record, strip_chr=False):
        # Add reference context
        if strip_chr:
            region = "{0}:{1}-{2}".format(
                vcf_record.chrom.replace("chr", "")
                if vcf_record.chrom != "chrM"
                else "MT",
                max(1, vcf_record.pos - self.context_size),
                vcf_record.stop + self.context_size,
            )
        else:
            region = "{0}:{1}-{2}".format(
                vcf_record.chrom,
                max(1, vcf_record.pos - self.context_size),
                vcf_record.stop + self.context_size,
            )
        context_record = get_builder(
            "CONTEXT", self.scheme, value=self.fa.fetch(region=region)
        )
        return context_record

    def shutdown(self):
        self.fa.close()


# __END__
