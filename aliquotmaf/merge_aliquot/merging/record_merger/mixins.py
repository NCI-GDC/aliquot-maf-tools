"""
Mixins for `BaseMafRecordMerger` classes
"""
from abc import ABCMeta, abstractmethod


class MafMergingAverageColumnsMixin(metaclass=ABCMeta):
    @abstractmethod
    def average_columns(self):
        """
        :return: a ``tuple`` of column names that should be averaged.
        """

    def do_mean_to_int(self, vals):
        """
        Rounds an array of values to the integer mean.

        :param vals: ``list`` of values to get average
        :return: the mean as an ``int``
        """
        return int("{0:.0f}".format(sum(vals) / float(len(vals))))


class MafMergingCombineColumnsMixin(metaclass=ABCMeta):
    @abstractmethod
    def combine_columns(self):
        """
        :return: a ``tuple`` of column names that should be combined into
        a unique set.
        """

    def do_uniq_list(self, records, column):
        """
        Gets a unique list of values.

        :param records: list of `maflib.MafRecord`s to extract data from
        :param column: the column key to merge
        :return: unique sorted list of values
        """
        vals = []
        for record in records:
            curr = record[column].value
            if isinstance(curr, list):
                vals.extend(curr)
            else:
                vals.append(curr)
        return sorted(list(set(vals)))
