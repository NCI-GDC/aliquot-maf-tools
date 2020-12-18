"""
Annotates the population frequency from the non-TCGA ExAC file.
"""
import pysam

from aliquotmaf.annotators.annotator import Annotator
from aliquotmaf.converters.builder import get_builder


class NonTcgaExac(Annotator):
    def __init__(self, scheme, source):
        super().__init__(name="NonTcgaExac", source=source, scheme=scheme)
        self.f = None
        self.popkeys = ["AFR", "AMR", "EAS", "FIN", "NFE", "OTH", "SAS"]

    @classmethod
    def setup(cls, scheme, source):
        curr = cls(scheme, source)
        curr.f = pysam.VariantFile(curr.source)
        return curr

    def annotate(self, maf_record, vcf_record, var_allele_idx=1):
        all_records = {}
        region = "{0}:{1}-{2}".format(vcf_record.chrom, vcf_record.pos, vcf_record.stop)
        alt = vcf_record.alleles[var_allele_idx]
        res = {}
        raw_af = None
        adj_af = None
        for record in self.f.fetch(region=region):
            if (
                vcf_record.pos == record.pos
                and vcf_record.ref == record.ref
                and alt in record.alts
            ):
                e_allele_idx = record.alts.index(alt)
                for p in self.popkeys:
                    ac_key = "AC_{0}".format(p)
                    an_key = "AN_{0}".format(p)
                    ac = record.info[ac_key][e_allele_idx]
                    an = record.info[an_key]
                    if an:
                        af = ac / float(an)
                        res[p] = af
                    else:
                        res[p] = None

                if record.info["AN"]:
                    raw_af = record.info["AC"][e_allele_idx] / float(record.info["AN"])
                if record.info["AN_Adj"]:
                    adj_af = record.info["AC_Adj"][e_allele_idx] / float(
                        record.info["AN_Adj"]
                    )
                break

        # Overall
        all_records["nontcga_ExAC_AF"] = get_builder(
            "nontcga_ExAC_AF", self.scheme, value=raw_af
        )
        all_records["nontcga_ExAC_AF_Adj"] = get_builder(
            "nontcga_ExAC_AF_Adj", self.scheme, value=adj_af
        )

        # pops
        for p in self.popkeys:
            key = "nontcga_ExAC_AF_{0}".format(p)
            all_records[key] = get_builder(key, self.scheme, value=res.get(p))

        return all_records

    def shutdown(self):
        self.f.close()
