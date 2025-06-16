"""
Builder classes assist in creating a MafColumnRecord.
"""

import abc
import sys

from maflib.column import MafColumnRecord


class Builder(metaclass=abc.ABCMeta):
    """
    Specify an abstract interface for creating column builder
    objects.
    """

    @classmethod
    @abc.abstractmethod
    def build(cls, column, value=None, default=None, scheme=None):
        """
        All builder classes must implement the build function.
        """


class BooleanColumnBuilder(Builder):
    """
    Builder class for BooleanColumn types.
    """

    @classmethod
    def build(cls, column, value=None, default=None, scheme=None):
        if value is None:
            assert default is not None, "BooleanColumn {0} can't be None".format(column)
            return MafColumnRecord.build(column, default, scheme=scheme)
        else:
            return MafColumnRecord.build(column, value, scheme=scheme)


class CanonicalBuilder(Builder):
    """
    Builder class for Canonical types.
    """

    @classmethod
    def build(cls, column, value=None, default=None, scheme=None):
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
    def build(cls, column, value=None, default=None, scheme=None):
        # None must be '' in this case
        default = default if default is not None else "None"
        if value is None:
            return MafColumnRecord.build(column, default, scheme=scheme)
        else:
            return MafColumnRecord.build(column, value, scheme=scheme)


class GenericColumnBuilder(Builder):
    """
    Generic builder class that should handle the majority of the cases.
    """

    @classmethod
    def build(cls, key, value, scheme, default=None, fn=None, **kwargs):
        if fn is not None:
            return MafColumnRecord.build(key, fn(value, **kwargs), scheme=scheme)
        elif value is None and default is not None:
            return MafColumnRecord.build(key, default, scheme=scheme)
        elif value is None and scheme.column_class(key).is_nullable():
            return scheme.column_class(key).build_nullable(key, scheme=scheme)
        else:
            return MafColumnRecord.build(key, value, scheme=scheme)


class GenericSequenceBuilder(Builder):
    """
    Generic sequence builder class that should handle the majority of the cases.
    """

    @classmethod
    def build(cls, key, value, scheme, default=None):
        if value is None and default is not None:
            return MafColumnRecord.build(key, default, scheme=scheme)
        elif value is None and scheme.column_class(key).is_nullable():
            return scheme.column_class(key).build_nullable(key, scheme=scheme)
        elif isinstance(value, list):
            return MafColumnRecord.build(key, ";".join(sorted(value)), scheme=scheme)
        else:
            return MafColumnRecord.build(key, value, scheme=scheme)


def get_builder(column, scheme, **kwargs):
    """
    Utility function to get the appropriate builder class.

    :param column: the column key
    :param scheme: the scheme class
    :returns: an appropriate builder class
    """
    colclassstr = scheme.column_class(column).__name__
    builderclassstr = "{0}Builder".format(colclassstr)
    try:
        if builderclassstr.startswith("SequenceOf"):
            builderclassstr = "GenericSequenceBuilder"
        return getattr(sys.modules[__name__], builderclassstr).build(
            column, scheme=scheme, **kwargs
        )
    except AttributeError:
        return GenericColumnBuilder.build(column, scheme=scheme, **kwargs)
