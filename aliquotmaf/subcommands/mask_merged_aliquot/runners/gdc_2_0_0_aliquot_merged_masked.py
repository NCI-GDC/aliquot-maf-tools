"""
Main logic for filtering merged aliquot MAFs for spec
gdc-2.0.0-aliquot-merged-masked.
"""
import json
from typing import Optional

from maflib.reader import MafReader
from maflib.record import MafRecord
from maflib.validation import ValidationStringency
from maflib.writer import MafWriter

from aliquotmaf.converters.utils import get_columns_from_header
from aliquotmaf.subcommands.mask_merged_aliquot.runners.gdc_1_0_0_aliquot_merged_masked import (
    GDC_1_0_0_Aliquot_Merged_Masked,
)

SPLICE_CONSEQUENCES = (
    'splice_acceptor_variant',
    'splice_donor_variant',
)


class GDC_2_0_0_Aliquot_Merged_Masked(GDC_1_0_0_Aliquot_Merged_Masked):
    def __init__(self, options: Optional[dict] = None):
        super(GDC_2_0_0_Aliquot_Merged_Masked, self).__init__(options)

        # Schema
        self.options["version"] = "gdc-1.0.0"
        self.options["annotation"] = "gdc-2.0.0-aliquot-merged-masked"

    def do_work(self) -> None:
        """Main wrapper function for running public MAF filter"""
        self.logger.info(
            "Processing input maf {0}...".format(self.options["input_maf"])
        )

        # Reader
        self.maf_reader = MafReader.reader_from(
            path=self.options["input_maf"],
            validation_stringency=ValidationStringency.Strict,
        )

        # Header
        self.setup_maf_header()

        # Writer
        self.maf_writer = MafWriter.from_path(
            path=self.options["output_maf"],
            header=self.maf_header,
            validation_stringency=ValidationStringency.Strict,
        )

        self._scheme = self.maf_header.scheme()  # type: ignore
        self._columns = get_columns_from_header(self.maf_header)
        self._colset = set(self._columns)

        # Counts
        processed = 0
        hotspot_gdc_set = set(["gdc_pon", "common_in_gnomAD"])
        nonexonic_set = set(["NonExonic"])

        record: MafRecord
        try:
            for record in self.maf_reader:  # type: ignore

                if processed > 0 and processed % 1000 == 0:
                    self.logger.info("Processed {0} records...".format(processed))

                callers = record["callers"].value  # type: ignore
                if (
                    len(callers) >= self.options["min_callers"]
                    and record["Mutation_Status"].value.value == "Somatic"  # type: ignore
                ):

                    self.metrics.add_sample_swap_metric(record)

                    gdc_filters = record["GDC_FILTER"].value  # type: ignore
                    gfset = set(gdc_filters)

                    if self.is_hotspot(record):
                        other_filts = gfset - hotspot_gdc_set
                        if len(other_filts) == 0:
                            self.write_record(record)
                        elif len(other_filts - nonexonic_set) == 0 and self.is_splice(
                            record
                        ):
                            # Rescue splicing if NonExonic
                            self.write_record(record)

                    # Rescue splicing if NonExonic
                    elif len(gfset - nonexonic_set) == 0 and self.is_splice(record):
                        self.write_record(record)

                    elif not gfset:
                        self.write_record(record)

                processed += 1
                self.metrics.input_records += 1

            self.logger.info("Processed {0} records.".format(processed))
            print(json.dumps(self.metrics.to_json(), indent=2, sort_keys=True))

        finally:

            self.maf_reader.close()
            self.maf_writer.close()

    def is_splice(self, record: MafRecord) -> bool:
        """
        Helper function that checks the One_Consequence column to check if the variant
        is a splice donor/acceptor.
        """
        cons = record["One_Consequence"].value  # type: ignore
        return cons in SPLICE_CONSEQUENCES

    @classmethod
    def __tool_name__(cls) -> str:
        return "gdc-2.0.0-aliquot-merged-masked"
