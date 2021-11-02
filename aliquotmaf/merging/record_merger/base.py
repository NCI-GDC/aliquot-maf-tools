"""
Base Class for merging MAF records.
"""
import logging
from abc import abstractmethod
from typing import TYPE_CHECKING, Any, Dict, List, Protocol

from aliquotmaf.converters.builder import get_builder
from aliquotmaf.logger import Logger

if TYPE_CHECKING:
    from maflib.column import MafColumnRecord
    from maflib.record import MafRecord
    from maflib.schemes import MafScheme

    from aliquotmaf.merging.overlap_set import OverlapSet


class BaseMafRecordMerger(Protocol):
    logger: logging.Logger
    scheme: 'MafScheme'
    columns: Any

    def __init__(self, scheme: 'MafScheme'):
        """
        Initialize the MAF record merging object which has the main `merge_records`
        function to take an `aliquotmaf.merging.overlap_set.OverlapSet`
        instance and performs merging.
        """

        self.logger = Logger.get_logger(self.__class__.__name__)
        self.scheme = scheme
        self.columns = scheme.column_names()

        self.logger.info("Loading MAF record merger...")

    @abstractmethod
    def caller_order(self) -> list:
        """
        :return: a ``list`` of caller names in their order of priority.
        """

    @abstractmethod
    def caller_type_order(self) -> list:
        """
        :return: a ``list`` of ``tuples`` of the format (caller, variant type)
        in their order of priority.
        """

    def allele_columns(self) -> tuple:
        """
        :return: a ``tuple`` of column names that contain allele information
        """
        return (
            "Tumor_Seq_Allele1",
            "Tumor_Seq_Allele2",
            "Match_Norm_Seq_Allele1",
            "Match_Norm_Seq_Allele2",
        )

    @abstractmethod
    def merge_records(
        self, results: 'OverlapSet', tumor_only: bool = False
    ) -> List['MafRecord']:
        """
        The main function for merging a MAF recrods from an
        `aliquotmaf.merging.overlap_set.OverlapSet` instance. The idea is to
        create a dictionary of column names as keys and the merged values as the
        values and then convert this dictionary to a `maflib.maf_record.MafRecord`
        instance and return it.

        :param results: an `aliquotmaf.merging.overlap_set.OverlapSet` instance
        :param tumor_only: ``True`` if there is no matched normal else ``False``
        :return: a list of merged `maflib.maf_record.MafRecord` instance
        """

    def standardize_alleles(
        self, maf_dic: Dict[str, 'MafColumnRecord'], tumor_only: bool = False
    ) -> Dict[str, 'MafColumnRecord']:
        """
        Helper utility to standardize all alleles to Ref/Alt tumor and
        Ref/Ref normal.

        :param maf_dic: ``dict`` of the maf record to format
        :param tumor_only: ``True`` if there is no matched normal else ``False``
        :return: updated maf_dic with formatted alleles
        """
        ref = maf_dic["Reference_Allele"].value
        alt = maf_dic["Allele"].value
        maf_dic["Tumor_Seq_Allele1"] = get_builder(
            "Tumor_Seq_Allele1", self.scheme, value=ref
        )
        maf_dic["Tumor_Seq_Allele2"] = get_builder(
            "Tumor_Seq_Allele2", self.scheme, value=alt
        )
        if tumor_only is False:
            norm_val = ref
        else:
            norm_val = None
        maf_dic["Match_Norm_Seq_Allele1"] = get_builder(
            "Match_Norm_Seq_Allele1", self.scheme, value=norm_val
        )
        maf_dic["Match_Norm_Seq_Allele2"] = get_builder(
            "Match_Norm_Seq_Allele2", self.scheme, value=norm_val
        )
        return maf_dic

    def fix_depths(
        self, maf_dic: Dict[str, 'MafColumnRecord'], tumor_only: bool = False
    ) -> Dict[str, 'MafColumnRecord']:
        """
        Sets the total depths of tumor/normal to be the sum of the ref and alt
        count columns if the sum of the ref and alt count columns is less than
        the dp.

        :param maf_dic: ``dict`` of the maf record to format
        :param tumor_only: ``True`` if there is no matched normal else ``False``
        :return: updated maf_dic with formatted depths
        """
        # fix depths
        tsum = maf_dic["t_ref_count"].value + maf_dic["t_alt_count"].value
        tdp = maf_dic["t_depth"].value
        if tsum > tdp:
            maf_dic["t_depth"] = get_builder("t_depth", self.scheme, value=tsum)

        if tumor_only is False:
            nsum = maf_dic["n_ref_count"].value + maf_dic["n_alt_count"].value
            ndp = maf_dic["n_depth"].value
            if nsum > ndp:
                maf_dic["n_depth"] = get_builder("n_depth", self.scheme, value=nsum)

        return maf_dic
