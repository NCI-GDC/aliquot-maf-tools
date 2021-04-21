"""
Main logic for filtering merged aliquot MAFs for spec
gdc-1.0.0-aliquot-merged-masked.
"""
import json

from maflib.header import MafHeader
from maflib.reader import MafReader
from maflib.sort_order import BarcodesAndCoordinate
from maflib.validation import ValidationStringency
from maflib.writer import MafWriter

from aliquotmaf.converters.builder import get_builder
from aliquotmaf.converters.utils import get_columns_from_header, init_empty_maf_record
from aliquotmaf.subcommands.mask_merged_aliquot.runners import BaseRunner


class GDC_1_0_0_Aliquot_Merged_Masked(BaseRunner):
    def __init__(self, options=dict()):
        super(GDC_1_0_0_Aliquot_Merged_Masked, self).__init__(options)

        # Schema
        self.options["version"] = "gdc-1.0.0"
        self.options["annotation"] = "gdc-1.0.0-aliquot-merged-masked"

    @classmethod
    def __add_arguments__(cls, parser):
        """Add the arguments to the parser"""
        parser.add_argument(
            "--tumor_only", action="store_true", help="Is this a tumor-only VCF?"
        )
        parser.add_argument(
            "--reference_fasta_index",
            required=False,
            help="Path to the reference fasta fai file if the input MAF is not sorted",
        )
        parser.add_argument(
            "--min_callers",
            default=2,
            type=int,
            help="Minimum number of callers required [2]",
        )

    def setup_maf_header(self):
        """
        Sets up the maf header.
        """
        # Reader header
        _hdr = MafHeader.from_reader(reader=self.maf_reader)

        if not self.options["reference_fasta_index"]:
            self.maf_header = MafHeader.from_defaults(
                version=self.options["version"],
                annotation=self.options["annotation"],
                sort_order=BarcodesAndCoordinate(),
                contigs=_hdr.contigs(),
            )
        else:
            self.maf_header = MafHeader.from_defaults(
                version=self.options["version"],
                annotation=self.options["annotation"],
                sort_order=BarcodesAndCoordinate(),
                fasta_index=self.options["reference_fasta_index"],
            )
        self.maf_header.validation_stringency = ValidationStringency.Strict

        header_date = BaseRunner.get_header_date()
        self.maf_header[header_date.key] = header_date

        try:
            nkey = _hdr["normal.aliquot"]
            self.maf_header["normal.aliquot"] = nkey
        except KeyError as e:
            if not self.options["tumor_only"]:
                raise e

        tkey = _hdr["tumor.aliquot"]
        self.maf_header["tumor.aliquot"] = tkey

    def do_work(self):
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

        self._scheme = self.maf_header.scheme()
        self._columns = get_columns_from_header(self.maf_header)
        self._colset = set(self._columns)

        # Counts
        processed = 0
        hotspot_gdc_set = set(["gdc_pon", "common_in_exac"])

        try:
            for record in self.maf_reader:

                if processed > 0 and processed % 1000 == 0:
                    self.logger.info("Processed {0} records...".format(processed))

                callers = record["callers"].value
                if (
                    len(callers) >= self.options["min_callers"]
                    and record["Mutation_Status"].value.value == "Somatic"
                ):

                    self.metrics.add_sample_swap_metric(record)

                    gdc_filters = record["GDC_FILTER"].value
                    gfset = set(gdc_filters)

                    if self.is_hotspot(record):
                        if len(gfset - hotspot_gdc_set) == 0:
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

    def is_hotspot(self, record):
        """
        Helper function to test if the record is marked as a hotspot.
        """
        if record["hotspot"].value and record["hotspot"].value.value == "Y":
            return True
        return False

    def write_record(self, record):
        """
        Helper function to write out the formatted merged public record.
        """
        self.metrics.collect_output(record)
        to_null = (
            "Match_Norm_Seq_Allele1",
            "Match_Norm_Seq_Allele2",
            "Match_Norm_Validation_Allele1",
            "Match_Norm_Validation_Allele2",
            "n_ref_count",
            "n_alt_count",
        )
        new_record = init_empty_maf_record()
        for column in self._columns:
            if column in to_null:
                new_record[column] = get_builder(column, self._scheme, value=None)
            else:
                new_record[column] = record[column]
        self.maf_writer += new_record

    @classmethod
    def __tool_name__(cls):
        return "gdc-1.0.0-aliquot-merged-masked"
