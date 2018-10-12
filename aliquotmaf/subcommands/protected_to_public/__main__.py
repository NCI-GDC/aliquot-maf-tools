"""
Subcommand for converting a protected aliquot MAF to a public aliquot MAF.
"""
from aliquotmaf.subcommands.base import Subcommand
from aliquotmaf.subcommands.protected_to_public.runners import GDC_1_2_0_Public

class ProtectedToPublic(Subcommand):
    @classmethod
    def __add_arguments__(cls, parser):
        """Add the arguments to the parser"""
        # Input group
        p_input = parser.add_argument_group(title="Input/Output Options")
        p_input.add_argument('--input_maf', required=True, help='Path to input protected MAF file')
        p_input.add_argument('--output_maf', required=True, help='Path to output public MAF file')
        
        # subparsers
        subparsers = parser.add_subparsers(dest="subcommand")
        subparsers.required = True

        GDC_1_2_0_Public.add(subparsers=subparsers)

    @classmethod
    def __get_description__(cls):
        """
        Optionally returns description
        """
        return "Filter and convert a protected aliquot-level MAF to a public aliquot-level MAF" 

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
