"""
Mixins for `BaseMafRecordMerger` classes
"""
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, List

if TYPE_CHECKING:
    from maflib.record import MafRecord


class MafMergingAverageColumnsMixin(ABC):
    @abstractmethod
    def average_columns(self) -> tuple:
        """
        :return: a ``tuple`` of column names that should be averaged.
        """

    def do_mean_to_int(self, vals: List[int]) -> int:
        """
        Rounds an array of values to the integer mean.

        :param vals: ``list`` of values to get average
        :return: the mean as an ``int``
        """
        return int(sum(vals) / len(vals))


class MafMergingCombineColumnsMixin(ABC):
    @abstractmethod
    def combine_columns(self) -> tuple:
        """
        :return: a ``tuple`` of column names that should be combined into
        a unique set.
        """

    def do_uniq_list(self, records: List['MafRecord'], column: str) -> list:
        """
        Gets a unique list of values.

        :param records: list of `maflib.MafRecord`s to extract data from
        :param column: the column key to merge
        :return: unique sorted list of values
        """
        vals = []
        for record in records:
            curr = record[column].value  # type: ignore
            if isinstance(curr, list):
                vals.extend(curr)
            else:
                vals.append(curr)
        return sorted(list(set(vals)))
