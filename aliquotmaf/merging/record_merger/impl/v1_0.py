"""
Maf recorder merger implementation v1.0
"""
from aliquotmaf.converters.builder import get_builder
from aliquotmaf.converters.utils import init_empty_maf_record
from aliquotmaf.merging.record_merger.base import BaseMafRecordMerger
from aliquotmaf.merging.record_merger.mixins import (
    MafMergingAverageColumnsMixin,
    MafMergingCombineColumnsMixin,
)


class MafRecordMerger_1_0_0(
    BaseMafRecordMerger, MafMergingAverageColumnsMixin, MafMergingCombineColumnsMixin
):
    def average_columns(self, tumor_only=False):
        """
        :return: a ``tuple`` of column names that should be averaged.
        """
        if not tumor_only:
            return (
                "t_depth",
                "t_ref_count",
                "t_alt_count",
                "n_depth",
                "n_ref_count",
                "n_alt_count",
            )
        else:
            return ("t_depth", "t_ref_count", "t_alt_count")

    def combine_columns(self):
        """
        :return: a ``tuple`` of column names that should be combined into
        a unique set.
        """
        return "GDC_FILTER"

    def caller_order(self):
        """
        :return: a ``list`` of caller names in their order of priority.
        """
        return [
            "vardict",
            "pindel",
            "mutect2",
            "muse",
            "varscan2",
            "caveman",
            "sanger_pindel",
            "gatk4_mutect2_pair",
            "somaticsniper",
        ]

    def caller_type_order(self):
        """
        :return: a ``list`` of ``tuples`` of the format (caller, variant type)
        in their order of priority.
        """
        # TODO: add svaba?
        return [
            ("mutect2", "MNP"),
            ("gatk4_mutect2_pair", "MNP"),
            ("vardict", "MNP"),
            ("pindel", "MNP"),
            ("sanger_pindel", "MNP"),
            ("caveman", "MNP"),
            ("mutect2", "DEL"),
            ("gatk4_mutect2_pair", "DEL"),
            ("vardict", "DEL"),
            ("pindel", "DEL"),
            ("sanger_pindel", "DEL"),
            ("varscan2", "DEL"),
            ("caveman", "DEL"),
            ("mutect2", "INS"),
            ("gatk4_mutect2_pair", "INS"),
            ("vardict", "INS"),
            ("pindel", "INS"),
            ("sanger_pindel", "INS"),
            ("varscan2", "INS"),
            ("caveman", "INS"),
            ("mutect2", "SNP"),
            ("gatk4_mutect2_pair", "SNP"),
            ("muse", "SNP"),
            ("vardict", "SNP"),
            ("varscan2", "SNP"),
            ("somaticsniper", "SNP"),
            ("caveman", "SNP"),
        ]

    def merge_records(self, results, tumor_only=False):
        """
        Returns one or more merged MAF records. This implementation is:

        1. If singleton -> format and write
        2. elif only one caller has overlaps -> loop over all records, format,
           and write
        3. elif there are multiple callers with only a single allele annotated
           -> merge into single record with most columns copied from caller
           based on self.caller_order()
        4. elif there are multiple variant types ->
               if there is a majority vote -> merge and write
               else -> collapse by self.caller_type_order (see self.collapse_by_caller_type)
        5. else -> collapse by self.caller_type_order (see self.collapse_by_caller_type)
        """
        maf_records = []

        if results.is_singleton():
            # Simple just extract first element from results for the caller and write
            maf_records = [
                self.maf_from_first_element(
                    results, results.callers, tumor_only=tumor_only
                )
            ]

        elif len(results.callers) == 1:
            # Simply will split the results for the single caller into separate maf records
            caller = results.callers[0]
            for variant in results[caller]:
                dic = {caller: [variant]}
                maf_records.append(
                    self.maf_from_first_element(
                        dic, results.callers, tumor_only=tumor_only
                    )
                )

        elif len(results.locus_allele_map) == 1:
            # Simply select the first element from all callers and create single record
            assert results.all_single_record()
            maf_records = [
                self.maf_from_first_element(
                    results, results.callers, tumor_only=tumor_only
                )
            ]

        else:
            if len(results.variant_types) == 1:
                curr = []
                _max = 1
                for allele in results.locus_allele_map:
                    ct = len(results.locus_allele_map[allele])
                    if ct > _max:
                        curr = [allele]
                        _max = ct
                    elif ct > 1 and ct == _max:
                        curr.append(allele)

                if len(curr) == 1:
                    val = results.locus_allele_map[curr[0]]
                    callers = sorted(list(val.keys()))
                    all_callers = set(results.callers)
                    star_callers = sorted(list(all_callers - set(callers)))
                    maf_records = [
                        self.maf_from_first_element(
                            val,
                            callers,
                            star_callers=star_callers,
                            tumor_only=tumor_only,
                        )
                    ]

                else:
                    maf_records = [
                        self.collapse_by_caller_type(results, tumor_only=tumor_only)
                    ]

            else:
                maf_records = [
                    self.collapse_by_caller_type(results, tumor_only=tumor_only)
                ]

        return maf_records

    def maf_from_first_element(
        self, results, callers, star_callers=[], tumor_only=False
    ):
        """
        Simply creates a MAF record from the first record in each caller.
        """
        maf_dic = {}
        selected_caller = None
        for caller in self.caller_order():
            if caller in callers:
                selected_caller = results[caller][0]
                break

        for column in self.columns:
            if column in self.allele_columns() or column == "callers":
                continue

            elif column in self.average_columns(tumor_only=tumor_only):
                vals = []
                for caller in callers:
                    curr = results[caller]
                    if not curr:
                        continue
                    curr = curr[0]
                    vals.append(curr[column].value)
                maf_dic[column] = get_builder(
                    column, self.scheme, value=self.do_mean_to_int(vals)
                )

            elif column in self.combine_columns():
                vals = self.do_uniq_list(
                    [results[i][0] for i in results if results[i]], column
                )
                maf_dic[column] = get_builder(column, self.scheme, value=vals)

            # NOTE: Not a solution, just a temp place holder until we fully build out the RNA annotator
            elif column == 'RNA_Support':
                maf_dic[column] = get_builder(column, self.scheme, value='Unknown')

            # NOTE: Not a solution, just a temp place holder until we fully build out the RNA annotator
            elif column in (
                'RNA_ref_count',
                'RNA_alt_count',
                'RNA_depth',
            ):
                maf_dic[column] = get_builder(column, self.scheme, value=None)

            else:
                maf_dic[column] = selected_caller[column]

        return self.format_dic_to_record(
            maf_dic, callers, star_callers=star_callers, tumor_only=tumor_only
        )

    def collapse_by_caller_type(self, results, tumor_only=False):
        """
        Use the caller type list to select.
        """
        all_callers = set(results.callers)
        selected_key = None

        for ctype in self.caller_type_order():
            if ctype in results.caller_type_map:
                selected_key = ctype
                break

        # If there is only a single record for the selected caller-type
        if len(results.caller_type_map[selected_key]) == 1:
            record = results.caller_type_map[selected_key][0]
            new_rec = {selected_key[0]: [record]}

            key = ":".join(
                list(
                    map(
                        str,
                        [
                            record["Start_Position"],
                            record["End_Position"],
                            record["Allele"],
                        ],
                    )
                )
            )
            # save lower priority matches
            other_matches = {
                k: results.locus_allele_map[key][k]
                for k in results.locus_allele_map[key]
                if k != selected_key[0]
            }
            new_rec.update(other_matches)
            star_callers = sorted(list(all_callers - set(new_rec.keys())))
            maf_record = self.maf_from_first_element(
                new_rec,
                sorted(list(new_rec.keys())),
                star_callers=star_callers,
                tumor_only=tumor_only,
            )
            return maf_record

        # When multiple overlapping records for selected caller-type, we
        # sort by vote and select first
        else:
            lst = []
            _selected_dic = {}
            for record in results.caller_type_map[selected_key]:
                key = ":".join(
                    list(
                        map(
                            str,
                            [
                                record["Start_Position"],
                                record["End_Position"],
                                record["Allele"],
                            ],
                        )
                    )
                )
                _selected_dic[key] = [record]
                n = len(results.locus_allele_map[key])
                val = (n, key)
                lst.append(val)
            selected_allele = sorted(lst, reverse=True)[0][1]
            new_rec = {selected_key[0]: _selected_dic[selected_allele]}
            other_matches = {
                k: results.locus_allele_map[selected_allele][k]
                for k in results.locus_allele_map[selected_allele]
                if k != selected_key[0]
            }
            new_rec.update(other_matches)
            star_callers = sorted(list(all_callers - set(new_rec.keys())))
            maf_record = self.maf_from_first_element(
                new_rec,
                sorted(list(new_rec.keys())),
                star_callers=star_callers,
                tumor_only=tumor_only,
            )
            return maf_record

    def format_dic_to_record(self, maf_dic, callers, star_callers=[], tumor_only=False):
        """
        Formats the dictionary into a MafRecord.
        """
        # Depths
        maf_dic = self.fix_depths(maf_dic, tumor_only=tumor_only)

        # Alleles
        maf_dic = self.standardize_alleles(maf_dic, tumor_only=tumor_only)

        # Callers
        _callers = callers + ["{0}*".format(i) for i in star_callers]
        maf_dic["callers"] = get_builder("callers", self.scheme, value=sorted(_callers))

        # Create MafRecord
        maf_record = init_empty_maf_record()
        for column in self.columns:
            idx = self.scheme.column_index(name=column)
            col = maf_dic[column]
            col.column_index = idx
            maf_record[column] = col
        return maf_record
