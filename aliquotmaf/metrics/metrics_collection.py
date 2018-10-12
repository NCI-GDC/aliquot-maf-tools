"""Class for collecting metrics from a MAF file."""
from collections import Counter

class MafMetricsCollection:
    def __init__(self):
        self.total_records = 0
        self.vcf_filters = Counter() 
        self.gdc_filters = Counter()
        self.variant_classification = Counter()
        self.variant_type = Counter()
        self.known = {'novel': 0, 'dbsnp': 0, 'cosmic': 0}

    def collect(self, record):
        """
        Collect metrics from the record.
        """
        self.total_records += 1
        # filters
        if record['FILTER'].value and record['FILTER'].value != ['PASS']:
            self.vcf_filters[str(record['FILTER'])] += 1

        if record['GDC_FILTER'].value: 
            self.gdc_filters[str(record['GDC_FILTER'])] += 1

        # variant clasification
        self.variant_classification[record['Variant_Classification'].value.value] += 1

        # variant type 
        self.variant_type[record['Variant_Type'].value.value] += 1

        # known
        if record['dbSNP_RS'].value:
            curr = record['dbSNP_RS'].value
            if curr == ['novel']:
                self.known['novel'] += 1
            else:
                self.known['dbsnp'] += 1

        if record['COSMIC'].value:
            self.known['cosmic'] += 1 

    def to_json(self):
        return {
            'total_records': self.total_records,
            'vcf_filters': self.vcf_filters,
            'gdc_filters': self.gdc_filters,
            'variant_classification': self.variant_classification,
            'variant_type': self.variant_type,
            'known': self.known
        }
