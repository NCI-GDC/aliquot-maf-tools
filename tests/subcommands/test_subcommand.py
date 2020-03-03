"""
Tests the subcommand ABC.
"""
import pytest

from aliquotmaf.subcommands.base import Subcommand
from aliquotmaf.__main__ import main


class Example(Subcommand):
    @classmethod
    def __add_arguments__(cls, subparser):
        pass

    @classmethod
    def __tool_name__(cls):
        return "example"


def test_get_title():
    assert Example.__tool_name__() == "example"


def test_get_description():
    assert Example.__get_description__() is None


def test_no_inputs():
    with pytest.raises(SystemExit):
        main(args=["example"])
