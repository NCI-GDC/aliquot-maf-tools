"""
Builder classes assist in creating a MafColumnRecord.
"""
import abc
import sys
from typing import TYPE_CHECKING, Any, Callable, Optional, Protocol

from maflib.column import MafColumnRecord

if TYPE_CHECKING:
    from maflib.schemes import MafScheme


# TODO: Refactor these classes. Should only change a default value
class Builder(Protocol):
    """
    Specify an abstract interface for creating column builder
    objects.
    """

    @classmethod
    @abc.abstractmethod
    def build(
        cls,
        column: str,
        value: Any = None,
        default: Any = None,
        scheme: Optional['MafScheme'] = None,
    ) -> MafColumnRecord:
        """
        All builder classes must implement the build function.
        """
        pass


class GenericColumnBuilder:
    """
    Generic builder class that should handle the majority of the cases.
    """

    @classmethod
    def build(
        cls,
        key: str,
        value: Any,
        scheme: 'MafScheme',
        default: Any = None,
        fn: Optional[Callable] = None,
        **kwargs: Any,
    ) -> MafColumnRecord:
        if fn is not None:
            return MafColumnRecord.build(key, fn(value, **kwargs), scheme=scheme)
        elif value is None and default is not None:
            return MafColumnRecord.build(key, default, scheme=scheme)
        elif value is None and scheme.column_class(key).is_nullable():
            return scheme.column_class(key).build_nullable(key, scheme=scheme)  # type: ignore
        else:
            return MafColumnRecord.build(key, value, scheme=scheme)


class BooleanColumnBuilder(Builder):
    """
    Builder class for BooleanColumn types.
    """

    @classmethod
    def build(
        cls,
        column: str,
        value: Any = None,
        default: Any = None,
        scheme: Optional['MafScheme'] = None,
    ) -> MafColumnRecord:
        if value is None:
            assert default is not None, f"BooleanColumn {column} can't be None"
            return MafColumnRecord.build(column, default, scheme=scheme)
        else:
            return MafColumnRecord.build(column, value, scheme=scheme)


class CanonicalBuilder(Builder):
    """
    Builder class for Canonical types.
    """

    @classmethod
    def build(
        cls,
        column: str,
        value: Any = None,
        default: Any = None,
        scheme: Optional['MafScheme'] = None,
    ) -> MafColumnRecord:
        # None must be '' in this case
        default = default if default is not None else ""
        if value is None:
            return MafColumnRecord.build(column, default, scheme=scheme)
        else:
            return MafColumnRecord.build(column, value, scheme=scheme)


class MutationStatusBuilder(Builder):
    """
    Builder class for MutationStatus types.
    """

    @classmethod
    def build(
        cls,
        column: str,
        value: Any = None,
        default: Any = None,
        scheme: Optional['MafScheme'] = None,
    ) -> MafColumnRecord:
        # None must be '' in this case
        default = default if default is not None else "None"
        if value is None:
            return MafColumnRecord.build(column, default, scheme=scheme)
        else:
            return MafColumnRecord.build(column, value, scheme=scheme)


class GenericSequenceBuilder(Builder):
    """
    Generic sequence builder class that should handle the majority of the cases.
    """

    @classmethod
    def build(
        cls,
        column: str,
        value: Any = None,
        default: Any = None,
        scheme: Optional['MafScheme'] = None,
    ) -> MafColumnRecord:
        if value is None and default is not None:
            return MafColumnRecord.build(column, default, scheme=scheme)
        elif (
            value is None
            and scheme is not None
            and scheme.column_class(column).is_nullable()
        ):
            return scheme.column_class(column).build_nullable(column, scheme=scheme)  # type: ignore
        elif isinstance(value, list):
            return MafColumnRecord.build(column, ";".join(sorted(value)), scheme=scheme)
        else:
            return MafColumnRecord.build(column, value, scheme=scheme)


def get_builder(column: str, scheme: 'MafScheme', **kwargs: Any) -> MafColumnRecord:
    """
    Utility function to get the appropriate builder class.

    :param column: the column key
    :param scheme: the scheme class
    :returns: an appropriate builder class
    """
    colclassstr = scheme.column_class(column).__name__
    builderclassstr = f"{colclassstr}Builder"
    try:
        if builderclassstr.startswith("SequenceOf"):
            builderclassstr = "GenericSequenceBuilder"
        return getattr(sys.modules[__name__], builderclassstr).build(  # type: ignore
            column, scheme=scheme, **kwargs
        )
    except AttributeError:
        return GenericColumnBuilder.build(column, scheme=scheme, **kwargs)
