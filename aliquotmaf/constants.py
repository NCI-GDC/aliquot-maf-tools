from typing import Final
from dataclasses import dataclass
import dataclasses

@dataclass()
class VariantCallerName:
    GDC_ENUM: str

    def snake(self):
        return self.GDC_ENUM.lower().replace(' ', '_')

    def option(self):
        return '--' + self.GDC_ENUM.lower().replace(' ', '-')

    def __repr__(self):
        return self.GDC_ENUM
    
    def __eq__(self, other):
        return self.__repr__() == other

@dataclass(frozen=True)
class VariantCallerConstants:
    MUTECT2: Final[VariantCallerName] = VariantCallerName("MuTect2")
    GATK4_MUTECT2: Final[VariantCallerName] = VariantCallerName("GATK4 MuTect2")
    GATK4_MUTECT2_PAIR: Final[VariantCallerName] = VariantCallerName("GATK4 MuTect2 Pair")
    SOMATIC_SNIPER: Final[VariantCallerName] = VariantCallerName("SomaticSniper")
    MUSE: Final[VariantCallerName] = VariantCallerName("MuSE")
    VARSCAN2: Final[VariantCallerName] = VariantCallerName("VarScan2")
    VARDICT: Final[VariantCallerName] = VariantCallerName("VarDict")
    CAVEMAN: Final[VariantCallerName] = VariantCallerName("CaVEMan")
    SANGER_PINDEL: Final[VariantCallerName] = VariantCallerName("Sanger Pindel")
    PINDEL: Final[VariantCallerName] = VariantCallerName("Pindel")
    SVABA_SOMATIC: Final[VariantCallerName] = VariantCallerName("SvABA Somatic")
    # provisional below
    STRELKA_SOMATIC: Final[VariantCallerName] = VariantCallerName("Strelka Somatic")

    def astuple(self):
        return tuple(flatten_tuple(dataclasses.astuple(self)))

def flatten_tuple(t):
    for x in t:
        if isinstance(x, tuple):
            yield from flatten_tuple(x)
        else:
            yield x 

variant_callers = VariantCallerConstants()

SPLICE_CONSEQUENCES = Final[set([
    'splice_acceptor_variant',
    'splice_donor_variant'
])]