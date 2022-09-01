"""A module for extracting allele information. All the classes present
here are subclasses of `~maf_converter_lib.extractor.Extractor` objects.

* ExtractVariantAlleleIndexParser   extracts the variant allele index
* ExtractGenotypeAndDepthsParser    extracts the genotype and depths
"""
from aliquotmaf.logger import Logger
from aliquotmaf.subcommands.vcf_to_aliquot.extractors import Extractor


class VariantAlleleIndexExtractor(Extractor):
    """Extractor class for extracting the variant allele index"""

    @classmethod
    def extract(cls, tumor_genotype):
        """
        Extracts the variant allele index from the tumor sample

        :param tumor_genotype: a dictionary or dictionary-like object
                               containing a GT key as a list or tuple
        :returns: the variant allele index
        """
        # First set var_allele_idx to 1
        var_allele_idx = 1
        # Get idx
        if tumor_genotype["GT"]:
            curr_gt = tumor_genotype["GT"]
            # Select the first non-REF allele
            try:
                var_allele_idx = [i for i in curr_gt if i is not None and i != 0][0]
            except IndexError:
                var_allele_idx = 1
        return var_allele_idx


class GenotypeAndDepthsExtractor(Extractor):
    """Extractor class for extracting the genotype and depths based on the
    variant allele index.
    """

    logger = Logger.get_logger("GenotypeAndDepthsExtractor")

    @classmethod
    def extract(cls, var_allele_idx, genotype, alleles):
        """
        Extracts the information for the variant alleles based on the
        variant allele index. Creates a new, updated genotype record
        and depths list.

        :param var_allele_idx: the variant allele index
        :param genotype: a dictionary or dictionary-like object containing
                         various possible keys like AD, DP, etc.
        :param alleles: an ordered list or tuple of the possible alleles
                        at the locus
        :returns: an updated genotype record and depths list
        """
        depths = []
        new_gt = {}
        if not genotype["GT"]:
            return new_gt, depths

        # If DP is defined, set it in new_gt
        if "DP" in genotype and genotype["DP"] is not None:
            new_gt["DP"] = genotype["DP"]

        # If AD is defined, then parse out all REF/ALT allele depths, or whatever is in it
        if "AD" in genotype and genotype["AD"] is not None:
            if isinstance(genotype["AD"], int):
                depths = [genotype["AD"]]
            else:
                depths = list(genotype["AD"])

        # handle VarScan VCF lines where AD contains only 1 depth, and REF allele depth is in RD
        if len(depths) == 1 and "RD" in genotype:
            depths = [None for i in alleles]
            depths[0] = genotype["RD"]
            if isinstance(genotype["AD"], int):
                depths[var_allele_idx] = genotype["AD"]
            else:
                depths[var_allele_idx] = genotype["AD"][0]

        # Handle SomaticSniper VCF lines, where allele depths must be extracted from BCOUNT
        elif "AD" not in genotype and "BCOUNT" in genotype:
            b_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
            bcount = list(genotype["BCOUNT"])
            depths = [bcount[b_idx[i]] if i in b_idx else None for i in alleles]

        # Handle CaVEMan which provides genotype and counts for all nucleotides in forward and reverse strands
        elif all(
            set(["GT", "FAZ", "FCZ", "FGZ", "FTZ", "RAZ", "RCZ", "RGZ", "RTZ"])
            in genotype.keys()
        ):
            var_allele = alleles[var_allele_idx]
            var_f_count_name = f"F{var_allele}Z"
            var_r_count_name = f"R{var_allele}Z"
            var_count = genotype[var_f_count_name] + genotype[var_r_count_name]
            ref_allele = [al for al in alleles if not var_allele][0]
            ref_f_count_name = f"F{ref_allele}Z"
            ref_r_count_name = f"R{ref_allele}Z"
            ref_count = genotype[ref_f_count_name] + genotype[ref_r_count_name]

            # set depths of alleles
            depths = [ref_count, var_count]

            # set DP to be sum of all observed base counts
            genotype["DP"] = sum(
                [
                    genotype[i]
                    for i in ["FAZ", "FCZ", "FGZ", "FTZ", "RAZ", "RCZ", "RGZ", "RTZ"]
                ]
            )

        # If N depths not equal to N alleles, blank out the depths
        elif depths and len(depths) != len(alleles):
            cls.logger.warning("The length of DP array != length of allele array")
            depths = [None for i in alleles]

        # Sanity check that REF/ALT allele depths are lower than total depth
        if (
            "DP" in genotype
            and genotype["DP"] is not None
            and (
                (depths[0] is not None and depths[0] > genotype["DP"])
                or (
                    depths[var_allele_idx] is not None
                    and depths[var_allele_idx] > genotype["DP"]
                )
                or (
                    depths[0] is not None
                    and depths[var_allele_idx] is not None
                    and depths[0] + depths[var_allele_idx] > genotype["DP"]
                )
            )
        ):
            cls.logger.warning("REF/ALT allele depths are lower than total depth!!")
            new_gt["DP"] = 0
            for i in depths:
                if i and i != ".":
                    new_gt["DP"] += i

        # If depths is empty, just set to 0, 0
        if not depths:
            depths = [0, 0]

        # If we have REF/ALT allele depths but not DP, then set DP equal to sum of all ADs
        if (depths[0] is not None and depths[var_allele_idx] is not None) and (
            "DP" not in genotype or genotype["DP"] is None or genotype["DP"] == "."
        ):
            # cls.logger.warn('Missing DP field. setting DP equal to sum of ADs!!')
            new_gt["DP"] = sum([i for i in depths if i and i != "."])

        # Set the formatted AD and alleles
        new_gt["AD"] = tuple([i if i != "" and i is not None else "." for i in depths])
        new_gt["GT"] = genotype["GT"]
        depths = [i if i != "." and i is not None else 0 for i in new_gt["AD"]]
        return new_gt, depths
