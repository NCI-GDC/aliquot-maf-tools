"""
Main entrypoint for all aliquot-maf-tools.
"""
import argparse
import datetime
import sys
from signal import SIG_DFL, SIGPIPE, signal

from aliquotmaf.logger import Logger
from aliquotmaf.subcommands.mask_merged_aliquot.__main__ import MaskMergedAliquotMaf
from aliquotmaf.subcommands.merge_aliquot.__main__ import MergeAliquotMafs
from aliquotmaf.subcommands.vcf_to_aliquot.__main__ import VcfToAliquotMaf

signal(SIGPIPE, SIG_DFL)


def main(args=None):
    """
    The main method for aliquot-maf-tools.
    """
    # Setup logger
    Logger.setup_root_logger()

    logger = Logger.get_logger("main")

    # Print header
    logger.info("-" * 75)
    logger.info("Program Args: aliquot-maf-tools " + " ".join(sys.argv[1::]))
    logger.info("Date/time: {0}".format(datetime.datetime.now()))
    logger.info("-" * 75)
    logger.info("-" * 75)

    # Get args
    p = argparse.ArgumentParser("GDC Aliquot MAF Tools")
    subparsers = p.add_subparsers(dest="subcommand")
    subparsers.required = True

    VcfToAliquotMaf.add(subparsers=subparsers)
    MergeAliquotMafs.add(subparsers=subparsers)
    MaskMergedAliquotMaf.add(subparsers=subparsers)

    options = p.parse_args(args)

    # Run
    cls = options.func(options)
    cls.do_work()

    # Finish
    logger.info("Finished!")


if __name__ == "__main__":
    main()
