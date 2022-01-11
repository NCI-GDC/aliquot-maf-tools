"""Class for collecting metrics from a MAF file."""
from collections import Counter
from typing import TYPE_CHECKING, TypedDict

if TYPE_CHECKING:
    from maflib.record import MafRecord


class MetricsDict(TypedDict):
    novel: int
    dbsnp: int
    cosmic: int
    common_in_exac: int


class SampleSwapsDict(TypedDict):
    total: int
    common_in_exac: int


class RecordDict(TypedDict):
    input_records: int
    output_records: int
    variant_classification: Counter
    variant_type: Counter
    known: MetricsDict
    sample_swaps: SampleSwapsDict


class MafMetricsCollection:
    def __init__(self) -> None:
        self.input_records: int = 0
        self.output_records: int = 0
        self.variant_classification: Counter = Counter()
        self.variant_type: Counter = Counter()
        self.known: MetricsDict = {
            "novel": 0,
            "dbsnp": 0,
            "cosmic": 0,
            "common_in_exac": 0,
        }
        self.sample_swaps: SampleSwapsDict = {"total": 0, "common_in_exac": 0}

    def add_sample_swap_metric(self, record: 'MafRecord') -> None:
        self.sample_swaps["total"] += 1

        # known
        is_common = "common_in_exac" in record["GDC_FILTER"].value  # type: ignore
        if is_common:
            self.sample_swaps["common_in_exac"] += 1

    def collect_output(self, record: 'MafRecord') -> None:
        """
        Collect metrics from the record that you will output.
        """
        self.output_records += 1

        # variant clasification
        self.variant_classification[record["Variant_Classification"].value.value] += 1  # type: ignore

        # variant type
        self.variant_type[record["Variant_Type"].value.value] += 1  # type: ignore

        # known
        is_common = "common_in_exac" in record["GDC_FILTER"].value  # type: ignore
        if is_common:
            self.known["common_in_exac"] += 1

        if record["dbSNP_RS"].value:  # type: ignore
            curr = record["dbSNP_RS"].value  # type: ignore
            if curr == ["novel"]:
                if not is_common:
                    self.known["novel"] += 1
            else:
                self.known["dbsnp"] += 1

        if record["COSMIC"].value:  # type: ignore
            self.known["cosmic"] += 1

    def to_json(self) -> RecordDict:
        return_dict: RecordDict = {
            "input_records": self.input_records,
            "output_records": self.output_records,
            "variant_classification": self.variant_classification,
            "variant_type": self.variant_type,
            "known": self.known,
            "sample_swaps": self.sample_swaps,
        }
        return return_dict
