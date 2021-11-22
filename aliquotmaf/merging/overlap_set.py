"""
Class containing a set of overlapping records and utilities.
"""

from typing import TYPE_CHECKING, Any, Dict, Iterator, List

if TYPE_CHECKING:
    from maflib.record import MafRecord


class OverlapSet:
    def __init__(self, result: List[List['MafRecord']], maf_keys: list):
        """
        Container for the results of an iteration of the overlap iterator. The
        results are converted to a dictionary where the caller is the key and the
        overlaps for that caller are the values.

        :param result: result from an iteration of `maflib.overlap_iter.LocatableOverlapIterator`
        :param maf_keys: list of callers in same order as results
        """
        self._data: Dict[str, Any] = dict(zip(maf_keys, result))
        self._has_indels = None

        self._callers: List[str]
        self._locus_allele_map: Dict[str, Dict[str, list]]
        self._caller_type_map: dict
        self._variant_types: tuple

        self._type_dic = {
            "SNP": "SNP",
            "DNP": "MNP",
            "TNP": "MNP",
            "ONP": "MNP",
            "DEL": "DEL",
            "INS": "INS",
        }

    def __iter__(self) -> Iterator[Any]:
        """
        Iterate over caller-records dictionary.
        """
        return iter(self._data)

    def __getitem__(self, key: str) -> Any:
        """
        Get results from particular caller like a dictionary.

        :param key: the caller id
        """
        return self._data[key]

    def is_singleton(self) -> bool:
        """
        Returns ``True`` if only a single variant from a single caller is
        present.
        """
        count = 0
        for caller in self._data:
            count += len(self._data[caller])
        if count == 1:
            return True
        return False

    @property
    def callers(self) -> List[str]:
        """
        Returns a list of callers containing an overlapping MAF record.
        """
        if getattr(self, "_callers", None) is None:
            lst = []
            for caller in sorted(self._data):
                if self._data[caller]:
                    lst.append(caller)
            self._callers = lst
        return self._callers

    @property
    def variant_types(self) -> tuple:
        """
        Returns a ``tuple`` of variant types detected. The variant types are mapped
        using the self._type_dic dictionary lookup which converts DNP/TNP/ONP to
        MNP.
        """
        if getattr(self, "_variant_types", None) is None:
            lst = []
            for caller in self:
                for record in self[caller]:
                    vtype = self._type_dic[record["Variant_Type"].value.value]
                    lst.append(vtype)
            self._variant_types = tuple(sorted(list(set(lst))))
        return self._variant_types

    @property
    def locus_allele_map(self) -> Dict[str, Dict[str, list]]:
        """
        Dictionary where the keys are Start_Postion:End_Position:Allele MAF columns
        and values are lists of records.
        """
        if getattr(self, "_locus_allele_map", None) is None:
            self._locus_allele_map = {}
            for caller in self._data:
                for record in self._data[caller]:
                    key = f"{record['Start_Position']}:{record['End_Position']}:{record['Allele']}"
                    if key not in self._locus_allele_map:
                        self._locus_allele_map[key] = {}
                    if caller not in self._locus_allele_map[key]:
                        self._locus_allele_map[key][caller] = []
                    self._locus_allele_map[key][caller].append(record)
        return self._locus_allele_map

    @property
    def caller_type_map(self) -> dict:
        """
        Dictionary where the keys are a ``tuple`` of caller and variant types
        mapped using self._type_dic and values are a list of records.
        """
        if getattr(self, "_caller_type_map", None) is None:
            self._caller_type_map = {}
            for caller in self._data:
                for record in self._data[caller]:
                    vtype = self._type_dic[record["Variant_Type"].value.value]
                    key = (caller, vtype)
                    if key not in self._caller_type_map:
                        self._caller_type_map[key] = []
                    self._caller_type_map[key].append(record)
        return self._caller_type_map

    def all_single_record(self) -> bool:
        """
        Returns ``True`` if all the callers contain one or zero records.
        """
        for caller in self.callers:
            if len(self[caller]) > 1:
                return False
        return True
