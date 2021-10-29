"""
Abstract base class for all subcommands in aliquot-maf-tool.
"""
from abc import ABCMeta, abstractmethod
from typing import TYPE_CHECKING, Protocol

if TYPE_CHECKING:
    from argparse import ArgumentParser, _SubParsersAction


class Subcommand(Protocol):
    @classmethod
    @abstractmethod
    def __add_arguments__(cls, parser: 'ArgumentParser') -> None:
        """Add the arguments to the parser"""

    @classmethod
    def __get_description__(cls) -> str:
        """
        Optionally returns description
        """

    @classmethod
    @abstractmethod
    def __tool_name__(cls) -> str:
        """
        Tool name to use for the subparser
        """

    @classmethod
    def add(cls, subparsers: '_SubParsersAction') -> 'ArgumentParser':
        """Adds the given subcommand to the subprsers."""
        subparser: 'ArgumentParser' = subparsers.add_parser(
            name=cls.__tool_name__(), description=cls.__get_description__()
        )

        cls.__add_arguments__(subparser)
        return subparser
