"""Utilities for converters"""
from typing import TYPE_CHECKING, List, Optional

from maflib.record import MafRecord
from maflib.validation import ValidationStringency

if TYPE_CHECKING:
    from maflib.header import MafHeader


def get_columns_from_header(header: 'MafHeader') -> List[str]:
    """
    Get the column names from the MAF header.
    """
    return header.scheme().column_names()  # type: ignore


def init_empty_maf_record(
    line_number: Optional[int] = None,
    stringency: ValidationStringency = ValidationStringency.Strict,
) -> MafRecord:
    """
    Initialize an empty maf record.
    """
    return MafRecord(line_number=line_number, validation_stringency=stringency)
