"""
MutationStatus annotator. Sets the Germline/Somatic/LOH/Unknown etc status for MAFs.
"""

from __future__ import absolute_import

from aliquotmaf.constants import variant_callers
from aliquotmaf.converters.builder import get_builder

from .annotator import Annotator


class MutationStatus(Annotator):
    def __init__(self, scheme, caller):
        super().__init__(name="MutationStatus", scheme=scheme)
        self.caller = caller
        self.mapper = {
            variant_callers.MUTECT2.GDC_ENUM: self._always_somatic,
            variant_callers.GATK4_MUTECT2.GDC_ENUM: self._always_somatic,
            variant_callers.SOMATIC_SNIPER.GDC_ENUM: self._somaticsniper,
            variant_callers.MUSE.GDC_ENUM: self._muse,
            variant_callers.VARSCAN2.GDC_ENUM: self._varscan,
            variant_callers.PINDEL.GDC_ENUM: self._muse,
            variant_callers.VARDICT.GDC_ENUM: self._vardict,
            variant_callers.CAVEMAN.GDC_ENUM: self._always_somatic,
            variant_callers.SANGER_PINDEL.GDC_ENUM: self._always_somatic,
            variant_callers.GATK4_MUTECT2_PAIR.GDC_ENUM: self._always_somatic,
            variant_callers.SVABA_SOMATIC.GDC_ENUM: self._always_somatic,
            variant_callers.STRELKA_SOMATIC.GDC_ENUM: self._always_somatic,  # needs validation
        }

    @classmethod
    def setup(cls, scheme, caller):
        curr = cls(scheme, caller)
        return curr

    def annotate(self, maf_record, vcf_record, tumor_sample):
        maf_record["Mutation_Status"] = get_builder(
            "Mutation_Status",
            self.scheme,
            value=self.mapper[self.caller](vcf_record, tumor_sample),
        )
        return maf_record

    def _always_somatic(self, record, tumor_sample):
        """Always somatic for MuTect2"""
        return "Somatic"

    def _somaticsniper(self, record, tumor_sample):
        """
        SomaticSniper parse from FORMAT and tumor GT column.
        """
        lookup = {0: "None", 1: "Germline", 2: "Somatic", 3: "LOH", 4: "Unknown"}
        curr = record.samples[tumor_sample]["SS"]
        return lookup[curr]

    def _varscan(self, record, tumor_sample):
        """VarScan2 parse from INFO column"""
        lookup = {
            "0": "None",
            "1": "Germline",
            "2": "Somatic",
            "3": "LOH",
            "5": "Unknown",
        }
        ss = record.info["SS"]
        return lookup[ss]

    def _muse(self, record, tumor_sample):
        """MuSE parser from FORMAT and tumor GT"""
        lookup = {
            0: "None",
            1: "Germline",
            2: "Somatic",
            3: "LOH",
            4: "Post-transcriptional modification",
            5: "Unknown",
        }
        ss = record.samples[tumor_sample]["SS"]
        return lookup[ss]

    def _vardict(self, record, tumor_sample):
        """VarDict parser from INFO"""
        lookup = {
            "LikelyLOH": "LOH",
            "StrongLOH": "LOH",
            "LikelySomatic": "Somatic",
            "Germline": "Germline",
            "AFDiff": "Unknown",
            "StrongSomatic": "Somatic",
        }
        val = record.info["STATUS"]
        return lookup[val]

    def shutdown(self):
        pass
