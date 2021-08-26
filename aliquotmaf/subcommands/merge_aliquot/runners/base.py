"""
Base class for all aliquot-level per-caller protected MAFs -> merged protected MAF runners.
"""
import datetime
from abc import ABCMeta, abstractmethod

from maflib.header import MafHeaderRecord

from aliquotmaf.logger import Logger


class BaseRunner(metaclass=ABCMeta):
    def __init__(self, options=dict()):
        self.logger = Logger.get_logger(self.__class__.__name__)
        self.options = options

        self.maf_readers = []
        self.callers = []
        self.maf_writer = None
        self._scheme = None
        self._columns = None
        self._colset = None
        self._merger = None

    @staticmethod
    def get_header_date():
        """
        Returns a MafHeaderRecord of the filedate.
        """
        return MafHeaderRecord(
            key="filedate", value=datetime.date.today().strftime("%Y%m%d")
        )

    @classmethod
    def __validate_options__(cls, options):
        """
        Optional function to validate other options
        """
        pass

    @classmethod
    def __get_description__(cls):
        """
        Optionally returns description
        """
        return None

    @classmethod
    def from_args(cls, args):
        cls.__validate_options__(args)
        return cls(options=vars(args))

    @abstractmethod
    def do_work(self):
        """Main wrapper function for running"""

    @classmethod
    @abstractmethod
    def __add_arguments__(cls, parser):
        """Add the arguments to the parser"""

    @classmethod
    @abstractmethod
    def __tool_name__(cls):
        """
        Tool name to use for the subparser
        """

    @classmethod
    def add(cls, subparsers):
        """Adds the given subcommand to the subparsers."""
        subparser = subparsers.add_parser(
            name=cls.__tool_name__(), description=cls.__get_description__()
        )

        cls.__add_arguments__(subparser)
        subparser.set_defaults(func=cls.from_args)
        return subparser
