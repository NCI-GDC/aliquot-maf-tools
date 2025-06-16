"""Utilities for converters"""

from maflib.record import MafRecord
from maflib.validation import ValidationStringency


def get_columns_from_header(header):
    """
    Get the column names from the MAF header.
    """
    return header.scheme().column_names()


def init_empty_maf_record(line_number=None, stringency=ValidationStringency.Strict):
    """
    Initialize an empty maf record.
    """
    return MafRecord(line_number=line_number, validation_stringency=stringency)
