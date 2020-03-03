"""
Subcommand for masking/filtering a merged aliquot MAF. 
"""
from aliquotmaf.subcommands.base import Subcommand
from aliquotmaf.subcommands.mask_merged_aliquot.runners import (
    GDC_1_0_0_Aliquot_Merged_Masked,
)


class MaskMergedAliquotMaf(Subcommand):
    @classmethod
    def __add_arguments__(cls, parser):
        """Add the arguments to the parser"""
        # Input group
        p_input = parser.add_argument_group(title="Input/Output Options")
        p_input.add_argument(
            "--input_maf", required=True, help="Path to input protected MAF file"
        )
        p_input.add_argument(
            "--output_maf", required=True, help="Path to output public MAF file"
        )

        # subparsers
        subparsers = parser.add_subparsers(dest="subcommand")
        subparsers.required = True

        GDC_1_0_0_Aliquot_Merged_Masked.add(subparsers=subparsers)

    @classmethod
    def __get_description__(cls):
        """
        Optionally returns description
        """
        return "Filter and mask a merged aliquot-level MAF"

    @classmethod
    def __tool_name__(cls):
        """
        Tool name to use for the subparser
        """
        return cls.__name__

    @classmethod
    def add(cls, subparsers):
        """Adds the given subcommand to the subparsers."""
        subparser = subparsers.add_parser(
            name=cls.__tool_name__(), description=cls.__get_description__()
        )

        cls.__add_arguments__(subparser)
        return subparser
