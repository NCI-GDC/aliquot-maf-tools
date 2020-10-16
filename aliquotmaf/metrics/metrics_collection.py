"""Class for collecting metrics from a MAF file."""
from collections import Counter


class MafMetricsCollection:
    def __init__(self):
        self.input_records = 0
        self.output_records = 0
        self.variant_classification = Counter()
        self.variant_type = Counter()
        self.known = {"novel": 0, "dbsnp": 0, "cosmic": 0, "common_in_exac": 0}
        self.sample_swaps = {"total": 0, "common_in_exac": 0}

    def add_sample_swap_metric(self, record):
        self.sample_swaps["total"] += 1

        # known
        is_common = "common_in_exac" in record["GDC_FILTER"].value
        if is_common:
            self.sample_swaps["common_in_exac"] += 1

    def collect_output(self, record):
        """
        Collect metrics from the record that you will output.
        """
        self.output_records += 1

        # variant clasification
        self.variant_classification[record["Variant_Classification"].value.value] += 1

        # variant type
        self.variant_type[record["Variant_Type"].value.value] += 1

        # known
        is_common = "common_in_exac" in record["GDC_FILTER"].value
        if is_common:
            self.known["common_in_exac"] += 1

        if record["dbSNP_RS"].value:
            curr = record["dbSNP_RS"].value
            if curr == ["novel"]:
                if not is_common:
                    self.known["novel"] += 1
            else:
                self.known["dbsnp"] += 1

        if record["COSMIC"].value:
            self.known["cosmic"] += 1

    def to_json(self):
        return {
            "input_records": self.input_records,
            "output_records": self.output_records,
            "variant_classification": self.variant_classification,
            "variant_type": self.variant_type,
            "known": self.known,
            "sample_swaps": self.sample_swaps,
        }
