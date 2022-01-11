"""
Base class for all aliquot-level per-caller protected MAFs -> merged protected MAF runners.
"""
import datetime
from abc import abstractmethod
from typing import TYPE_CHECKING, Any, List, Optional, Protocol

from maflib.header import MafHeaderRecord

from aliquotmaf.logger import Logger

if TYPE_CHECKING:
    from argparse import ArgumentParser, Namespace, _SubParsersAction

    from maflib.reader import MafReader
    from maflib.schemes import MafScheme
    from maflib.writer import MafWriter


class BaseRunner(Protocol):
    logger: Logger
    options: dict
    callers: list
    maf_readers: List['MafReader']
    maf_writer: 'MafWriter'
    _scheme: 'MafScheme'
    _columns: List[str]
    _colset: set
    _merger: Any

    def __init__(self, options: Optional[dict] = None):
        self.logger = Logger.get_logger(self.__class__.__name__)
        self.options = options if options is not None else {}

        self.maf_readers: List['MafReader'] = []

    @staticmethod
    def get_header_date() -> MafHeaderRecord:
        """
        Returns a MafHeaderRecord of the filedate.
        """
        return MafHeaderRecord(
            key="filedate", value=datetime.date.today().strftime("%Y%m%d")
        )

    @abstractmethod
    def do_work(self) -> None:
        """Main wrapper function for running"""

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
    def __get_description__(cls) -> Optional[str]:
        """
        Optionally returns description
        """
        return None

    @classmethod
    def __validate_options__(cls, options: 'Namespace') -> None:
        pass

    @classmethod
    def from_args(cls, args: 'Namespace') -> 'BaseRunner':
        cls.__validate_options__(args)
        return cls(options=vars(args))

    @classmethod
    def add(cls, subparsers: '_SubParsersAction') -> 'ArgumentParser':
        """Adds the given subcommand to the subparsers."""
        subparser = subparsers.add_parser(
            name=cls.__tool_name__(), description=cls.__get_description__()
        )

        cls.__add_arguments__(subparser)
        subparser.set_defaults(func=cls.from_args)
        return subparser
