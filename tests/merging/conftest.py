from collections import OrderedDict

import pytest
from maflib.column_types import (
    DnaString,
    NullableDnaString,
    NullableStringColumn,
    NullableZeroBasedIntegerColumn,
    OneBasedIntegerColumn,
    SequenceOfStrings,
    StringColumn,
    VariantType,
    YesNoOrUnknown,
    ZeroBasedIntegerColumn,
)
from maflib.record import MafRecord
from maflib.validation import ValidationStringency

from aliquotmaf.merging.overlap_set import OverlapSet


@pytest.fixture
def test_output_scheme(get_test_scheme):
    vals = [
        ("Chromosome", StringColumn),
        ("Start_Position", OneBasedIntegerColumn),
        ("End_Position", OneBasedIntegerColumn),
        ("Variant_Type", VariantType),
        ("Reference_Allele", DnaString),
        ("Allele", DnaString),
        ("Tumor_Seq_Allele1", DnaString),
        ("Tumor_Seq_Allele2", DnaString),
        ("Match_Norm_Seq_Allele1", NullableDnaString),
        ("Match_Norm_Seq_Allele2", NullableDnaString),
        ("t_depth", ZeroBasedIntegerColumn),
        ("t_ref_count", ZeroBasedIntegerColumn),
        ("t_alt_count", ZeroBasedIntegerColumn),
        ("n_depth", NullableZeroBasedIntegerColumn),
        ("n_ref_count", NullableZeroBasedIntegerColumn),
        ("n_alt_count", NullableZeroBasedIntegerColumn),
        ("extra", NullableStringColumn),
        ("GDC_FILTER", SequenceOfStrings),
        ("RNA_Support", YesNoOrUnknown),
        ("RNA_depth", NullableZeroBasedIntegerColumn),
        ("RNA_ref_count", NullableZeroBasedIntegerColumn),
        ("RNA_alt_count", NullableZeroBasedIntegerColumn),
        ("callers", SequenceOfStrings),
    ]

    coldict = OrderedDict(vals)
    return get_test_scheme(
        coldict, version="test-merging-out", annotation="test-merging-out-anno"
    )


@pytest.fixture
def test_input_scheme(get_test_scheme):
    vals = [
        ("Chromosome", StringColumn),
        ("Start_Position", OneBasedIntegerColumn),
        ("End_Position", OneBasedIntegerColumn),
        ("Variant_Type", VariantType),
        ("Reference_Allele", DnaString),
        ("Allele", DnaString),
        ("Tumor_Seq_Allele1", DnaString),
        ("Tumor_Seq_Allele2", DnaString),
        ("Match_Norm_Seq_Allele1", NullableDnaString),
        ("Match_Norm_Seq_Allele2", NullableDnaString),
        ("t_depth", ZeroBasedIntegerColumn),
        ("t_ref_count", ZeroBasedIntegerColumn),
        ("t_alt_count", ZeroBasedIntegerColumn),
        ("n_depth", NullableZeroBasedIntegerColumn),
        ("n_ref_count", NullableZeroBasedIntegerColumn),
        ("n_alt_count", NullableZeroBasedIntegerColumn),
        ("extra", NullableStringColumn),
        ("GDC_FILTER", SequenceOfStrings),
    ]

    coldict = OrderedDict(vals)
    return get_test_scheme(coldict)


@pytest.fixture
def overlapped_records_generate():
    def _generate_overlaps(test_input_scheme, line_list, callers):
        lst = []

        for line in line_list:
            if isinstance(line, list):
                curr = []
                for item in line:
                    curr.append(
                        MafRecord.from_line(
                            item,
                            scheme=test_input_scheme,
                            validation_stringency=ValidationStringency.Strict,
                        )
                    )
                lst.append(curr)
            elif line:
                lst.append(
                    [
                        MafRecord.from_line(
                            line,
                            scheme=test_input_scheme,
                            validation_stringency=ValidationStringency.Strict,
                        )
                    ]
                )
            else:
                lst.append([])

        return OverlapSet(lst, callers)

    return _generate_overlaps
