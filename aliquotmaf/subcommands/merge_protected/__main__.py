"""
Subcommand for merging single caller protected MAFs to merged protected MAFs. 
"""
from aliquotmaf.subcommands.base import Subcommand
from aliquotmaf.subcommands.merge_protected.runners import GDC_1_0_0_Merged_Protected

class MergeProtected(Subcommand):
    @classmethod
    def __add_arguments__(cls, parser):
        """Add the arguments to the parser"""
        # Input group
        p_input = parser.add_argument_group(title="Input/Output Options")
        p_input.add_argument('--output_maf', required=True, help='Path to output public MAF file')
        
        # subparsers
        subparsers = parser.add_subparsers(dest="subcommand")
        subparsers.required = True

        GDC_1_0_0_Merged_Protected.add(subparsers=subparsers)

    @classmethod
    def __get_description__(cls):
        """
        Optionally returns description
        """
        return "Merge protected MAFs from multiple callers within the same aliquot"

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
            name=cls.__tool_name__(),
            description=cls.__get_description__())

        cls.__add_arguments__(subparser)
        return subparser
