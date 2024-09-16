from typing import Final
from dataclasses import dataclass
import dataclasses

@dataclass(frozen=True)
class VariantCallerNamespace:
    MUTECT2: Final[str] = "MuTect2"
    GATK4_MUTECT2: Final[str] = "GATK4 MuTect2"
    SOMATIC_SNIPER: Final[str] = "SomaticSniper"
    MUSE: Final[str] = "MuSE"
    VARSCAN2: Final[str] = "VarScan2"
    PINDEL: Final[str] = "Pindel"
    VARDICT: Final[str] = "VarDict"
    CAVEMAN: Final[str] = "CaVEMan"
    SANGER_PINDEL: Final[str] = "Sanger Pindel"
    GATK4_MUTECT2_PAIR: Final[str] = "GATK4 MuTect2 Pair"
    SVABA: Final[str] = "SvABA"
    # provisional below
    STRELKA2_MANTA: Final[str] = "Strelka2 Manta"

    def astuple(cls):
        return dataclasses.astuple(cls)

variant_callers = VariantCallerNamespace()