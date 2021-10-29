"""
A simple container class for storing a collection of builders you want to convert to
MAF columns.
"""
from enum import Enum
from typing import TYPE_CHECKING, Any, Generator, List, Optional

from maflib.schemes import MafScheme

from aliquotmaf.converters.builder import get_builder

if TYPE_CHECKING:
    from maflib.column import MafColumnRecord
    from maflib.schemes import MafScheme


class StateEnum(Enum):
    TRANSFORMED = "TRANSFORMED"
    UNTRANSFORMED = "UNTRANSFORMED"


class InputData:
    """Contains individual column mapping and transformation functions."""

    def __init__(
        self, value: Optional[str], column: str, default: Optional[str] = None
    ):
        self.value = value
        self.column = column
        self.default = default
        self._transformed: 'MafColumnRecord'
        self.state: str = StateEnum.UNTRANSFORMED.value

    # FIXME: Overwriting builtin, bad
    def __hash__(self) -> str:  # type: ignore
        return self.column

    def __str__(self) -> str:
        return f"{self.state}: <FROM({self.value})> -> <TO({self.column} -> {self.transformed})>"

    def __repr__(self) -> str:
        return str(self)

    def __build__(self, scheme: 'MafScheme') -> None:
        self._transformed = get_builder(
            self.column, scheme, value=self.value, default=self.default
        )
        self.state = StateEnum.TRANSFORMED.value

    @property
    def transformed(self) -> Optional['MafColumnRecord']:
        if getattr(self, "_transformed", None) is not None:
            return self._transformed
        return None


class InputCollection:
    """
    Holds the data you want to transform into MafColumnRecords. You add builders
    to the instance and then use transform to convert to MAF columns.
    """

    def __init__(self) -> None:
        self._from_data: List[InputData] = []
        self._colset: dict = {}

    def columns(self) -> list:
        return [i.column for i in self._from_data]

    def add(
        self,
        value: Optional[str],
        column: str,
        default: Optional[str] = None,
        data: Optional[InputData] = None,
    ) -> None:
        if data is not None:
            assert (
                data.column not in self._colset
            ), f"Column '{column}' is already present!"
            self._colset[data.column] = None
            self._from_data.append(data)
        else:
            assert column not in self._colset, f"Column '{column}' is already present!"
            self._colset[column] = None
            self._from_data.append(InputData(value, column, default=default))

    def transform(self, scheme: 'MafScheme') -> None:
        for i in self._from_data:
            i.__build__(scheme)

    def __iter__(self) -> Generator[Any, None, None]:
        for i in self._from_data:
            yield i


# __END__
