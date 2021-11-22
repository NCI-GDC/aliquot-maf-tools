"""
Base class for all protected -> public MAF runners.
"""
import datetime
import logging
from abc import abstractmethod
from typing import TYPE_CHECKING, Optional, Protocol

from maflib.header import MafHeaderRecord

from aliquotmaf.logger import Logger
from aliquotmaf.metrics.metrics_collection import MafMetricsCollection

if TYPE_CHECKING:
    from argparse import ArgumentParser, Namespace, _SubParsersAction

    from maflib.header import MafHeader
    from maflib.reader import MafReader
    from maflib.schemes import MafScheme
    from maflib.writer import MafWriter


class BaseRunner(Protocol):

    logger: logging.Logger
    options: dict
    maf_reader: 'MafReader'
    maf_writer: 'MafWriter'
    maf_header: 'MafHeader'
    metrics: MafMetricsCollection
    _scheme: 'MafScheme'
    _columns: list
    _colset: set

    @staticmethod
    def get_header_date() -> MafHeaderRecord:
        """
        Returns a MafHeaderRecord of the filedate.
        """
        return MafHeaderRecord(
            key="filedate", value=datetime.date.today().strftime("%Y%m%d")
        )

    def __init__(self, options: Optional[dict] = None):
        self.logger = Logger.get_logger(self.__class__.__name__)
        self.options = options if options is not None else {}

        self.metrics = MafMetricsCollection()

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
        pass

    @classmethod
    def __get_description__(cls) -> Optional[str]:
        """
        Optionally returns description
        """
        pass

    @classmethod
    def from_args(cls, args: 'Namespace') -> 'BaseRunner':
        cls.__validate_options__(args)
        return cls(options=vars(args))

    @classmethod
    def __validate_options__(cls, options: 'Namespace') -> None:
        pass

    @classmethod
    def add(cls, subparsers: '_SubParsersAction') -> 'ArgumentParser':
        """Adds the given subcommand to the subparsers."""
        subparser = subparsers.add_parser(
            name=cls.__tool_name__(), description=cls.__get_description__()
        )

        cls.__add_arguments__(subparser)
        subparser.set_defaults(func=cls.from_args)
        return subparser
