"""A module for extracting allele information. All the classes present
here are subclasses of `~maf_converter_lib.extractor.Extractor` objects.

* ExtractVariantAlleleIndexParser   extracts the variant allele index
* ExtractGenotypeAndDepthsParser    extracts the genotype and depths
"""
from typing import Final, List

from aliquotmaf.logger import Logger
from aliquotmaf.subcommands.vcf_to_aliquot.extractors import Extractor

# Variant Caller Enumeration
# Theoretical possible values:
# - CaVEMan
# - GATK4 MuTect2
# - MuSE
# - MuTect2
# - Pindel
# - Sanger Pindel
# - SomaticSniper
# - VarDict
# - VarScan2
# - Strelka2 RNA
# - SvABA Somatic
# - GATK4 MuTect2 Pair

# Values for which specific extraction algorithms have been implemented
SANGER_CAVEMAN: Final[str] = "CaVEMan"
SANGER_PINDEL: Final[str] = "Sanger Pindel"

CAVEMAN_NUCLEOTIDE_COUNTS: Final[List[str]] = [
    "FAZ",
    "FCZ",
    "FGZ",
    "FTZ",
    "RAZ",
    "RCZ",
    "RGZ",
    "RTZ",
]

PINDEL_NUCLEOTIDE_COUNTS: Final[List[str]] = [
    "ND",
    "PD",
    "NR",
    "PR",
    "NB",
    "PB",
    "NP",
    "PP",
    "NU",
    "PU",
    "FD",
    "FC",
]


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

    # This should be re-written to separate the logic used to extract these data
    # from different variant caller outputs
    @classmethod
    def extract(cls, var_allele_idx, genotype, alleles, caller_id):
        """
        Extracts the information for the variant alleles based on the
        variant allele index. Creates a new, updated genotype record
        and depths list.

        :param var_allele_idx: the variant allele index
        :param genotype: a dictionary or dictionary-like object containing
                         various possible keys like AD, DP, etc.
        :param alleles: an ordered list or tuple of the possible alleles
                        at the locus
        :returns: an updated genotype record and formatted depths list
        """
        depths = []
        new_gt = {}
        if caller_id == SANGER_PINDEL:
            new_gt, depths = cls.extract_sanger_pindel(
                var_allele_idx, genotype, alleles
            )
        elif caller_id == SANGER_CAVEMAN:
            new_gt, depths = cls.extract_sanger_caveman(
                var_allele_idx, genotype, alleles
            )
        else:
            new_gt, depths = cls.extract_legacy(var_allele_idx, genotype, alleles)
            if new_gt == {} and depths == []:
                return new_gt, depths
        ad, f_depths = cls.format_results(depths)
        new_gt["AD"] = ad
        return new_gt, f_depths

    @classmethod
    def format_results(cls, depths):
        """
        Format outputs for AD and depths
        """
        # Set the formatted AD and alleles
        ad = tuple([i if i != "" and i is not None else "." for i in depths])
        nd = [i if i != "." and i is not None else 0 for i in ad]
        return ad, nd

    @classmethod
    def extract_legacy(cls, var_allele_idx, genotype, alleles):
        depths = []
        new_gt = {}
        # TODO: Add test for this branch
        if not genotype["GT"]:
            return new_gt, depths

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

        # # Handle CaVEMan which provides genotype and counts for all nucleotides in forward and reverse strands
        # elif set(["GT"] + CAVEMAN_NUCLEOTIDE_COUNTS).issubset(genotype.keys()):
        #     new_gt, depths = cls.extract_sanger_caveman(var_allele_idx, genotype, alleles)

        # # Handle Sanger Pindel
        # elif set(["GT"] + PINDEL_NUCLEOTIDE_COUNTS).issubset(genotype.keys()):
        #     new_gt, depths = cls.extract_sanger_pindel(var_allele_idx, genotype, alleles)

        # If N depths not equal to N alleles, blank out the depths
        elif depths and len(depths) != len(alleles):
            cls.logger.warning("The length of DP array != length of allele array")
            depths = [None for i in alleles]

        # If DP is defined, set it in new_gt
        if "DP" in genotype and genotype["DP"] is not None:
            new_gt["DP"] = genotype["DP"]

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
            cls.logger.warning("REF/ALT allele depths are larger than total depth!!")
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
        new_gt["GT"] = genotype["GT"]
        return new_gt, depths

    @classmethod
    def extract_sanger_caveman(cls, var_allele_idx, genotype, alleles):
        """
        Extracts depths, counts, and genotype information for sanger pindel vcfs

        :param var_allele_idx: the variant allele index
        :param genotype: a dictionary or dictionary-like object containing
                         various possible keys like AD, DP, etc.
        :param alleles: an ordered list or tuple of the possible alleles
                        at the locus
        :returns: an updated genotype record and depths list

        Sanger CaVEMan output contains counts for all possible nucleotides
        in forward and reverse strands.

        in order to extract counts we have to first get the ref and alt nucleotides
        and then if for example the desired ref or alt nucleotide is 'G' construct
        a key like 'FGZ' and 'RGZ' to access the forward and reverse counts of G.
        """

        var_allele = alleles[var_allele_idx]
        var_f_count_name = f"F{var_allele}Z"
        var_r_count_name = f"R{var_allele}Z"
        var_count = genotype[var_f_count_name] + genotype[var_r_count_name]
        ref_allele = [al for al in alleles if al != var_allele][0]
        ref_f_count_name = f"F{ref_allele}Z"
        ref_r_count_name = f"R{ref_allele}Z"
        ref_count = genotype[ref_f_count_name] + genotype[ref_r_count_name]

        # set depths of alleles
        depths = [ref_count, var_count]

        # set DP to be sum of all observed base counts
        new_gt = {}
        new_gt["DP"] = sum([genotype[i] for i in CAVEMAN_NUCLEOTIDE_COUNTS])
        new_gt["GT"] = genotype["GT"]
        return (new_gt, depths)

    @classmethod
    def extract_sanger_pindel(cls, var_allele_idx, genotype, alleles):
        """
        Extracts depths, counts, and genotype information for sanger pindel vcfs

        :param var_allele_idx: the variant allele index
        :param genotype: a dictionary or dictionary-like object containing
                         various possible keys like AD, DP, etc.
        :param alleles: an ordered list or tuple of the possible alleles
                        at the locus
        :returns: an updated genotype record and depths list
        """

        depths = []
        new_gt = {}
        total_depth = genotype['NR'] + genotype['PR']
        alt_count = genotype['NP'] + genotype['PP']
        ref_count = total_depth - alt_count
        new_gt["DP"] = total_depth
        new_gt["GT"] = genotype["GT"]
        depths = [ref_count, alt_count]
        return (new_gt, depths)
