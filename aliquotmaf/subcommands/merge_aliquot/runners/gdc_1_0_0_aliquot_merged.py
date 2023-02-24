"""
Main logic for merging raw aliquot MAFs on schema gdc-1.0.0-aliquot-merged.
"""

from maflib.header import MafHeader
from maflib.overlap_iter import LocatableOverlapIterator
from maflib.reader import MafReader
from maflib.sort_order import BarcodesAndCoordinate
from maflib.sorter import MafSorter
from maflib.validation import ValidationStringency
from maflib.writer import MafWriter

import aliquotmaf.filters as Filters
from aliquotmaf.converters.builder import get_builder
from aliquotmaf.converters.utils import get_columns_from_header
from aliquotmaf.merging.filtering_iterator import FilteringPeekableIterator
from aliquotmaf.merging.overlap_set import OverlapSet
from aliquotmaf.merging.record_merger.impl.v1_0 import MafRecordMerger_1_0_0
from aliquotmaf.subcommands.merge_aliquot.runners import BaseRunner


class GDC_1_0_0_Aliquot_Merged(BaseRunner):
    def __init__(self, options=dict()):
        super(GDC_1_0_0_Aliquot_Merged, self).__init__(options)

        # Schema
        self.options["version"] = "gdc-1.0.0"
        self.options["annotation"] = "gdc-1.0.0-aliquot-merged"

    @classmethod
    def __add_arguments__(cls, parser):
        """Add the arguments to the parser"""
        parser.add_argument(
            "--tumor_only", action="store_true", help="If this is a tumor-only MAF"
        )
        parser.add_argument(
            "--mutect2", help="Path to input protected MuTect2 MAF file"
        )
        parser.add_argument("--muse", help="Path to input protected MuSE MAF file")
        parser.add_argument(
            "--vardict", help="Path to input protected VarDict MAF file"
        )
        parser.add_argument(
            "--varscan2", help="Path to input protected VarScan2 MAF file"
        )
        parser.add_argument(
            "--somaticsniper", help="Path to input protected SomaticSniper MAF file"
        )
        parser.add_argument("--pindel", help="Path to input protected Pindel MAF file")
        parser.add_argument(
            "--min_n_depth",
            type=int,
            default=7,
            help="Flag variants where normal depth is <= INT as ndp. "
            + "This is performed after averaging "
            + "depths across callers [7]",
        )
        parser.add_argument(
            "--caveman", help="Path to input protected CaVEMan MAF file"
        )
        parser.add_argument(
            "--sanger-pindel", help="Path to input protected Sanger Pindel MAF file"
        )
        parser.add_argument(
            "--gatk4-mutect2-pair",
            help="Path to input protected GATK4 MuTect2 Pair MAF file",
        )
        parser.add_argument(
            "--gatk4-mutect2", help="Path to input protected GATK4 MuTect2 MAF file"
        )

    def load_readers(self):
        """
        Loads the array of MafReaders and sets the callers list.
        """
        # TODO: Add more callers
        maf_keys = [
            "mutect2",
            "muse",
            "vardict",
            "varscan2",
            "somaticsniper",
            "pindel",
            "caveman",
            "sanger_pindel",
            "gatk4_mutect2_pair",
            "gatk4_mutect2",
        ]

        for maf_key in maf_keys:
            if self.options[maf_key]:
                self.logger.info("{0} MAF {1}".format(maf_key, self.options[maf_key]))
                self.maf_readers.append(
                    MafReader.reader_from(
                        path=self.options[maf_key],
                        validation_stringency=ValidationStringency.Strict,
                    )
                )
                self.callers.append(maf_key)

    def setup_maf_header(self):
        """
        Sets up the maf header.
        """
        # Reader header
        _hdr = MafHeader.from_reader(reader=self.maf_readers[0])

        self.maf_header = MafHeader.from_defaults(
            version=self.options["version"],
            annotation=self.options["annotation"],
            sort_order=BarcodesAndCoordinate(),
            contigs=_hdr.contigs(),
        )
        self.maf_header.validation_stringency = ValidationStringency.Strict

        header_date = BaseRunner.get_header_date()
        self.maf_header[header_date.key] = header_date

        if not self.options["tumor_only"]:
            nkey = _hdr["normal.aliquot"]
            self.maf_header["normal.aliquot"] = nkey
        tkey = _hdr["tumor.aliquot"]
        self.maf_header["tumor.aliquot"] = tkey

    def do_work(self):
        """Main wrapper function for running protect MAF merging"""

        # Reader
        self.load_readers()

        # Header
        self.setup_maf_header()

        self._scheme = self.maf_header.scheme()
        self._columns = get_columns_from_header(self.maf_header)

        # Sorter
        sorter = MafSorter(
            max_objects_in_ram=100000,
            sort_order_name=BarcodesAndCoordinate.name(),
            scheme=self.maf_header.scheme(),
            contigs=self.maf_header.contigs(),
        )

        # Merger
        self._merger = MafRecordMerger_1_0_0(self._scheme)

        # Overlap iterator
        o_iter = LocatableOverlapIterator(
            self.maf_readers,
            contigs=self.maf_header.contigs(),
            peekable_iterator_class=FilteringPeekableIterator,
        )

        # Set up Normal Depth Filter for tumor-only or tumor-normal operation
        if self.options['tumor_only']:
            ndp_filter = Filters.NormalDepth.setup(None)
        else:
            ndp_filter = Filters.NormalDepth.setup(self.options["min_n_depth"])
        ndp_tag = ndp_filter.tags[0]

        # Counts
        processed = 0
        try:
            for record in o_iter:
                # progress update
                if processed > 0 and processed % 1000 == 0:
                    self.logger.info(
                        "Processed {0} overlapping intervals...".format(processed)
                    )

                result = OverlapSet(record, self.callers)

                for maf_record in self._merger.merge_records(
                    result, tumor_only=self.options["tumor_only"]
                ):
                    if maf_record is not None:
                        # Recheck normal depth
                        gdc_filters = maf_record["GDC_FILTER"].value
                        has_tag = ndp_tag in gdc_filters
                        ndp = ndp_filter.filter(maf_record)
                        if has_tag != ndp:
                            if ndp:
                                gdc_filters.extend(ndp_filter.tags)
                            else:
                                gdc_filters = list(
                                    filter(
                                        lambda x: x != ndp_filter.tags[0], gdc_filters
                                    )
                                )

                            maf_record["GDC_FILTER"] = get_builder(
                                "GDC_FILTER", self._scheme, value=sorted(gdc_filters)
                            )

                        # Add to sorter
                        sorter += maf_record

                processed += 1

            self.logger.info("Writing {0} sorted, merged records...".format(processed))

            # Writer
            self.maf_writer = MafWriter.from_path(
                path=self.options["output_maf"],
                header=self.maf_header,
                validation_stringency=ValidationStringency.Strict,
            )

            counter = 0
            for record in sorter:
                if counter > 0 and counter % 1000 == 0:
                    self.logger.info(
                        "Wrote {0} sorted, merged records...".format(counter)
                    )
                self.maf_writer += record
                counter += 1

            self.logger.info(
                "Finished writing {0} sorted, merged records.".format(counter)
            )

        finally:
            for reader in self.maf_readers:
                reader.close()

            sorter.close()

            if self.maf_writer:
                self.maf_writer.close()

    @classmethod
    def __tool_name__(cls):
        return "gdc-1.0.0-aliquot-merged"
