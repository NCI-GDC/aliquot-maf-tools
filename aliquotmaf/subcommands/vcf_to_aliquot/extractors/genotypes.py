"""A module for extracting allele information. All the classes present
here are subclasses of `~maf_converter_lib.extractor.Extractor` objects.

* ExtractVariantAlleleIndexParser   extracts the variant allele index
* ExtractGenotypeAndDepthsParser    extracts the genotype and depths
"""
from typing import Final, List

from aliquotmaf.constants import variant_callers
from aliquotmaf.logger import Logger
from aliquotmaf.subcommands.vcf_to_aliquot.extractors import Extractor

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

    @classmethod
    def extract(cls, var_allele_idx, genotype, alleles, caller_id):
        """
        Decides which extraction algorithm to use based on the caller_id.

        :param var_allele_idx: the variant allele index
        :param genotype: a dictionary or dictionary-like object containing
                         various possible keys like AD, DP, etc.
        :param alleles: an ordered list or tuple of the possible alleles
                        at the locus
        :param caller_id: a string identifying the variant caller that 
                          generated the vcf file
        :returns: an updated genotype record and depths list
        """
        depths = []
        new_gt = {}
        if not genotype["GT"]:
            return new_gt, depths
        
        if caller_id == variant_callers.CAVEMAN:
            depths, dp = cls.__extract_caveman(var_allele_idx, genotype, alleles)
        if caller_id == variant_callers.GATK4_MUTECT2:
            depths, dp = cls.__extract_mutect2(genotype)
        if caller_id == variant_callers.MUTECT2:
            depths, dp = cls.__extract_mutect2(genotype)
        if caller_id == variant_callers.GATK4_MUTECT2_PAIR:
            depths, dp = cls.__extract_mutect2(genotype)
        if caller_id == variant_callers.SOMATIC_SNIPER:
            depths, dp = cls.__extract_somaticsniper(var_allele_idx, genotype, alleles)
        if caller_id == variant_callers.MUSE:
            depths, dp = cls.__extract_muse(genotype)
        if caller_id == variant_callers.VARSCAN2:
            depths, dp = cls.__extract_varscan2(var_allele_idx, genotype, alleles)
        if caller_id == variant_callers.PINDEL:
            depths, dp = cls.__extract_pindel(genotype)
        if caller_id == variant_callers.SANGER_PINDEL:
            depths, dp = cls.__extract_sanger_pindel(genotype)
        if caller_id == variant_callers.SVABA_SOMATIC:
            depths, dp = cls.__extract_svaba_somatic(genotype)
        if caller_id == variant_callers.STRELKA_SOMATIC:
            depths, dp = cls.__strelka_somatic(genotype)
        
        # Set the formatted AD and alleles
        new_gt["AD"] = tuple([i if i != "" and i is not None else "." for i in depths])
        new_gt["GT"] = genotype["GT"]
        new_gt["DP"] = dp
        depths = [i if i != "." and i is not None else 0 for i in new_gt["AD"]]
        return new_gt, depths



            # if caller_id == variant_callers.VARDICT:
            #     self.extract_legacy(cls, var_allele_idx, genotype, alleles)


            # if caller_id == variant_callers.STRELKA2_MANTA:
            #     self.extract_legacy(cls, var_allele_idx, genotype, alleles)

            # if caller_id == variant_callers.ABSOLUTE:
            #     self.extract_legacy(cls, var_allele_idx, genotype, alleles)

    @classmethod
    def __extract_mutect2(cls, genotype):
        '''
        MuTect2 should always have the AD tag set and may or may not have 
        the DP flag set. This function handles those cases, adding up 
        allele-depths when necessary.
        '''

        DP_IS_SET = ("DP" in genotype and genotype["DP"] is not None and genotype["DP"] != '.')
        # set allele depths
        if isinstance(genotype["AD"], int):
            allele_depths = [genotype["AD"]]
        else:
            allele_depths = list(genotype["AD"])
        # set total depth
        if DP_IS_SET:
            dp = genotype["DP"]
        else:
            dp = sum([i for i in allele_depths if i and i != "."])
        return (allele_depths, dp)
    
    @classmethod
    def __extract_muse(cls, genotype):
        '''
        MuSE reports AD and DP tags
        '''
        if isinstance(genotype["AD"], int):
            allele_depths = [genotype["AD"]]
        else:
            allele_depths = list(genotype["AD"])
        dp = genotype["DP"]
        return (allele_depths, dp)

    @classmethod
    def __extract_caveman(cls, var_allele_idx, genotype, alleles):
        '''
        CaVEMan provides genotype and counts for all nucleotides
        in forward and reverse strand so we need to add them up
        to get allele depths and overall depth
        '''
        # Handle CaVEMan which provides genotype and counts for all nucleotides in forward and reverse strands
        var_allele = alleles[var_allele_idx]
        var_f_count_name = f"F{var_allele}Z"
        var_r_count_name = f"R{var_allele}Z"
        var_count = genotype[var_f_count_name] + genotype[var_r_count_name]
        ref_allele = [al for al in alleles if al != var_allele][0]
        ref_f_count_name = f"F{ref_allele}Z"
        ref_r_count_name = f"R{ref_allele}Z"
        ref_count = genotype[ref_f_count_name] + genotype[ref_r_count_name]

        # set depths of alleles
        allele_depths = [ref_count, var_count]

        # set DP to be sum of all observed base counts
        dp = sum([genotype[i] for i in CAVEMAN_NUCLEOTIDE_COUNTS])
        return (allele_depths, dp)

    @classmethod
    def __extract_somaticsniper(cls, var_allele_idx, genotype, alleles):
        '''
        SomaticSniper reports total depth as DP but allele depths must be
        extracted from BCOUNT
        '''
        # set allele-depths
        b_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
        bcount = list(genotype["BCOUNT"])
        allele_depths = [bcount[b_idx[i]] if i in b_idx else None for i in alleles]

        DP_IS_SET = ("DP" in genotype and genotype["DP"] is not None and genotype["DP"] != '.')
        # set total depth
        if DP_IS_SET:
            dp = genotype["DP"]
        else:
            dp = sum([i for i in allele_depths if i and i != "."])
        return (allele_depths, dp)

    @classmethod
    def __extract_varscan2(cls, var_allele_idx, genotype, alleles):
        '''
        VarScan2 only reports alt allele depth in 'AD'
        reference allele depth is in 'RD' tag
        overall depth is in 'DP'
        '''
                # handle VarScan VCF lines where AD contains only 1 depth, and REF allele depth is in RD

        allele_depths = [None for i in alleles]
        allele_depths[0] = genotype["RD"]
        if isinstance(genotype["AD"], int):
            allele_depths[var_allele_idx] = genotype["AD"]
        else:
            allele_depths[var_allele_idx] = genotype["AD"][0]
        dp = genotype["DP"]
        return (allele_depths, dp)
    
    @classmethod
    def __extract_pindel(cls, genotype):
        '''
        Pindel reports AD and DP tags
        '''
        if isinstance(genotype["AD"], int):
            allele_depths = [genotype["AD"]]
        else:
            allele_depths = list(genotype["AD"])
        dp = genotype["DP"]
        return (allele_depths, dp)
    
    @classmethod
    def __extract_sanger_pindel(cls, genotype):
        """
        Sanger pindel reports an assortment of counts on reverse and 
        forward strands. For Alt allele we use sum of 'Pindel' calls on 
        both strands. Total depth is take from sum of 'Total' mapped reads
        on both strands, and Ref allele depths is the difference of these.

        :param genotype: a dictionary or dictionary-like object containing
                         various possible keys like AD, DP, etc.

        :returns: an updated genotype record and depths list
        """
        total_depth = genotype['NR'] + genotype['PR']
        alt_count = genotype['NP'] + genotype['PP']
        ref_count = total_depth - alt_count
        dp = total_depth
        allele_depths = [ref_count, alt_count]
 
        return (allele_depths, dp)

    @classmethod
    def __extract_svaba_somatic(cls, genotype):
        """
        SvABA only reports the depth of the alt allele in 'AD' and reports
        total depth in 'DP' so we infer the ref allele depth from these.
        """
        total_depth = genotype['DP']
        alt_depth = genotype['AD']
        ref_depth = total_depth - alt_depth
        allele_depths = [ref_depth, alt_depth]
        dp = total_depth
        return (allele_depths, dp)

    def _extract_strelka_somatic(cls, genotype):
        """
        There are two subtypes we have to deal with 

        Strelka reports total depth in DP
        """
        pass

    def extract_legacy(cls, var_allele_idx, genotype, alleles):
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
        elif set(["GT"] + CAVEMAN_NUCLEOTIDE_COUNTS).issubset(genotype.keys()):
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
            new_gt["DP"] = sum([genotype[i] for i in CAVEMAN_NUCLEOTIDE_COUNTS])

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
                ([0] is not None and depths[0] > genotype["DP"])
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

        # Set the formatted AD and alleles
        new_gt["AD"] = tuple([i if i != "" and i is not None else "." for i in depths])
        new_gt["GT"] = genotype["GT"]
        depths = [i if i != "." and i is not None else 0 for i in new_gt["AD"]]
        return new_gt, depths
