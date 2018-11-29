"""
Class containing a set of overlapping records and utilities.
"""

class OverlapSet:
    def __init__(self, result, maf_keys):
        self._data = dict(zip(maf_keys, result))
        self._locus_allele_map = None
        self._caller_type_map = None
        self._has_indels = None
        self._callers = None
        self._variant_types = None

        self._type_dic = {
            'SNP': 'SNP', 'DNP': 'MNP', 'TNP': 'MNP', 'ONP': 'MNP', 'DEL': 'DEL', 'INS': 'INS'
        }

    def __iter__(self):
        return iter(self._data)

    def __getitem__(self, key):
        return self._data[key]

    def is_singleton(self):
        count = 0
        for caller in self._data:
            count += len(self._data[caller])
        if count == 1: return True
        return False

    @property
    def callers(self):
        if self._callers is None:
            lst = []
            for caller in sorted(self._data):
                if self._data[caller]:
                    lst.append(caller)
            self._callers = lst
        return self._callers

    @property
    def variant_types(self):
        if self._variant_types is None:
            lst = []
            for caller in self:
                for record in self[caller]:
                    vtype = self._type_dic[record['Variant_Type'].value.value]
                    lst.append(vtype)
            self._variant_types = tuple(sorted(list(set(lst))))
        return self._variant_types

    @property
    def locus_allele_map(self):
        if self._locus_allele_map is None:
            self._locus_allele_map = {}
            for caller in self._data:
                for record in self._data[caller]:
                    key = ':'.join(list(map(str, [record['Start_Position'], record['End_Position'],
                        record['Allele']])))
                    if key not in self._locus_allele_map: self._locus_allele_map[key] = {}
                    if caller not in self._locus_allele_map[key]: self._locus_allele_map[key][caller] = []
                    self._locus_allele_map[key][caller].append(record)
        return self._locus_allele_map

    @property
    def caller_type_map(self):
        if self._caller_type_map is None:
            self._caller_type_map = {}
            for caller in self._data:
                for record in self._data[caller]:
                    vtype = self._type_dic[record['Variant_Type'].value.value]
                    key = (caller, vtype)
                    if key not in self._caller_type_map:
                        self._caller_type_map[key] = []
                    self._caller_type_map[key].append(record)
        return self._caller_type_map

    def all_single_record(self):
        for caller in self.callers:
            if len(self[caller]) > 1:
                return False
        return True
