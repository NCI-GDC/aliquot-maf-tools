"""
MutationStatus annotator. Sets the Germline/Somatic/LOH/Unknown etc status for MAFs.
"""
from __future__ import absolute_import

import pysam

from .annotator import Annotator

from aliquotmaf.converters.builder import get_builder

class MutationStatus(Annotator):
    def __init__(self, scheme, caller):
        super().__init__(name='MutationStatus', scheme=scheme)
        self.caller = caller
        self.mapper = {
            'MuTect2': self._mutect2,
            'GATK4 MuTect2': self._mutect2,
            'SomaticSniper': self._somaticsniper,
            'MuSE': self._muse,
            'VarScan2': self._varscan
        }

    @classmethod
    def setup(cls, scheme, caller):
        curr = cls(scheme, caller)
        return curr

    def annotate(self, maf_record, vcf_record, tumor_sample):
        maf_record['Mutation_Status'] = get_builder(
            "Mutation_Status", 
            self.scheme, 
            value=self.mapper[self.caller](vcf_record, tumor_sample))
        return maf_record 


    def _mutect2(self, record, tumor_sample):
        """Always somatic for MuTect2"""
        return 'Somatic'

    def _somaticsniper(self, record, tumor_sample):
        """
        SomaticSniper parse from FORMAT and tumor GT column.
        """
        lookup = {0: 'None', 1: 'Germline', 2: 'Somatic', 3: 'LOH', 4: 'Unknown'}
        curr = record.samples[tumor_sample]['SS']
        return lookup[curr]

    def _varscan(self, record, tumor_sample):
        """VarScan2 parse from INFO column"""
        lookup = {'0': 'None', '1': 'Germline', '2': 'Somatic', '3': 'LOH', '5': 'Unknown'}
        ss = record.info['SS']
        return lookup[ss]

    def _muse(self, record, tumor_sample):
        """MuSE parser from FORMAT and tumor GT"""
        lookup = {0: 'None', 1: 'Germline', 2: 'Somatic', 
            3: 'LOH', 4: 'Post-transcriptional modification', 5: 'Unknown'}
        ss = record.samples[tumor_sample]['SS']
        return lookup[ss]

    def shutdown(self):
        pass 
