"""
Subcommand for converting a VEP annotated VCF to a protected
aliquot MAF.
"""
from aliquotmaf.subcommands.base import Subcommand
from aliquotmaf.subcommands.vcf_to_protected.runners import GDC_1_2_0_Protected

class VcfToProtected(Subcommand):
    @classmethod
    def __add_arguments__(cls, parser):
        """Add the arguments to the parser"""
        # Input group
        p_input = parser.add_argument_group(title="Input/Output Options")
        p_input.add_argument('--input_vcf', required=True, help='Path to input VCF file')
        p_input.add_argument('--output_maf', required=True, help='Path to output MAF file')
        
        # subparsers
        subparsers = parser.add_subparsers(dest="subcommand")
        subparsers.reqqqquired = True

        GDC_1_2_0_Protected.add(subparsers=subparsers)

    @classmethod
    def __get_description__(cls):
        """
        Optionally returns description
        """
        return "Convert VEP annotated VCF to an aliquot-level protected MAF" 

    @classmethod
    def __tool_name__(cls):
        """
        Tool name to use for the subparser
        """
        return cls.__name__

    @classmethod
    def add(cls, subparsers):
        """Adds the given subcommand to the subprsers."""
        subparser = subparsers.add_parser(
            name=cls.__tool_name__(),
            description=cls.__get_description__())

        cls.__add_arguments__(subparser)
        return subparser
