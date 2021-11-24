"""
Subcommand for converting a VEP annotated VCF to a raw aliquot MAF.
"""
from typing import TYPE_CHECKING

from aliquotmaf.subcommands.base import Subcommand
from aliquotmaf.subcommands.vcf_to_aliquot.runners.gdc_1_0_0_aliquot import (
    GDC_1_0_0_Aliquot,
)
from aliquotmaf.subcommands.vcf_to_aliquot.runners.gdc_2_0_0_aliquot import (
    GDC_2_0_0_Aliquot,
)

if TYPE_CHECKING:
    from argparse import ArgumentParser, _SubParsersAction


class VcfToAliquotMaf(Subcommand):
    @classmethod
    def __add_arguments__(cls, parser: 'ArgumentParser') -> None:
        """Add the arguments to the parser"""
        # Input group
        p_input = parser.add_argument_group(title="Input/Output Options")
        p_input.add_argument(
            "--input_vcf", required=True, help="Path to input VCF file"
        )
        p_input.add_argument(
            "--output_maf", required=True, help="Path to output MAF file"
        )

        # subparsers
        subparsers = parser.add_subparsers(dest="subcommand")
        subparsers.required = True

        GDC_1_0_0_Aliquot.add(subparsers=subparsers)
        GDC_2_0_0_Aliquot.add(subparsers=subparsers)

    @classmethod
    def __get_description__(cls) -> str:
        """
        Optionally returns description
        """
        return "Convert VEP annotated VCF to an aliquot-level raw MAF"

    @classmethod
    def __tool_name__(cls) -> str:
        """
        Tool name to use for the subparser
        """
        return cls.__name__

    @classmethod
    def add(cls, subparsers: '_SubParsersAction') -> 'ArgumentParser':
        """Adds the given subcommand to the subprsers."""
        subparser = subparsers.add_parser(
            name=cls.__tool_name__(), description=cls.__get_description__()
        )

        cls.__add_arguments__(subparser)
        return subparser
