"""
A simple container class for storing a collection of builders you want to convert to
MAF columns.
"""

from aliquotmaf.converters.builder import get_builder


class InputCollection:
    """
    Holds the data you want to transform into MafColumnRecords. You add builders
    to the instance and then use transform to convert to MAF columns.
    """

    def __init__(self):
        self._from_data = []
        self._colset = {}

    def columns(self):
        return [i.column for i in self._from_data]

    def add(self, value=None, column=None, default=None, data=None):
        if data is not None:
            assert data.column not in self._colset, (
                "Column '{0}' is already present!".format(data.column)
            )
            self._colset[data.column] = None
            self._from_data.append(data)
        else:
            assert column not in self._colset, (
                "Column '{0}' is already present!".format(column)
            )
            self._colset[column] = None
            self._from_data.append(InputData(value, column, default=default))

    def transform(self, scheme):
        for i in self._from_data:
            i.__build__(scheme)

    def __iter__(self):
        for i in self._from_data:
            yield i


class InputData:
    """Contains individual column mapping and transformation functions."""

    def __init__(self, value, column, default=None):
        self.value = value
        self.column = column
        self.default = default
        self.transformed = None
        self.state = "UNTRANSFORMED"

    def __hash__(self):
        return self.column

    def __str__(self):
        return "{0.state}: <FROM({0.value})> -> <TO({0.column} -> {0.transformed})>".format(
            self
        )

    def __repr__(self):
        return str(self)

    def __build__(self, scheme):
        self.transformed = get_builder(
            self.column, scheme, value=self.value, default=self.default
        )
        self.state = "TRANSFORMED"
