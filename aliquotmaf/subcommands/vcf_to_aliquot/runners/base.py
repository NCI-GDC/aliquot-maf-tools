"""Base class for all vcf2maf runners"""
import datetime
from abc import ABCMeta, abstractmethod
from typing import TYPE_CHECKING, List, Optional, Protocol

from maflib.header import MafHeaderRecord
from maflib.schemes import MafScheme

from aliquotmaf.logger import Logger

if TYPE_CHECKING:
    from argparse import ArgumentParser, Namespace, _SubParsersAction

    from maflib.header import MafHeader
    from maflib.schemes import MafScheme
    from maflib.writer import MafWriter


class BaseRunner(Protocol):
    logger: Logger
    options: dict
    maf_header: 'MafHeader'
    maf_writer: 'MafWriter'
    _scheme: 'MafScheme'
    _columns: List[str]
    _colset: set

    def __init__(self, options: Optional[dict] = None):
        self.logger = Logger.get_logger(self.__class__.__name__)
        self.options = options if options is not None else {}

        # Maf stuff
        self._scheme: 'MafScheme'
        self._columns: List[str]
        self._colset: set

    @staticmethod
    def get_header_date() -> MafHeaderRecord:
        """
        Returns a MafHeaderRecord of the filedate.
        """
        return MafHeaderRecord(
            key="filedate", value=datetime.date.today().strftime("%Y%m%d")
        )

    @classmethod
    def __validate_options__(cls, options: 'Namespace') -> None:
        """
        Optional function to validate other options
        """
        pass

    @classmethod
    def __get_description__(cls) -> Optional[str]:
        """
        Optionally returns description
        """
        return None

    @classmethod
    def from_args(cls, args: 'Namespace') -> 'BaseRunner':
        cls.__validate_options__(args)
        return cls(options=vars(args))

    @abstractmethod
    def do_work(self) -> None:
        """Main wrapper function for running vcf2maf"""

    @classmethod
    @abstractmethod
    def __add_arguments__(cls, parser: 'ArgumentParser') -> None:
        """Add the arguments to the parser"""

    @classmethod
    @abstractmethod
    def __tool_name__(cls) -> str:
        """
        Tool name to use for the subparser
        """

    @classmethod
    def add(cls, subparsers: '_SubParsersAction') -> 'ArgumentParser':
        """Adds the given subcommand to the subprsers."""
        subparser = subparsers.add_parser(
            name=cls.__tool_name__(), description=cls.__get_description__()
        )

        cls.__add_arguments__(subparser)
        subparser.set_defaults(func=cls.from_args)
        return subparser
