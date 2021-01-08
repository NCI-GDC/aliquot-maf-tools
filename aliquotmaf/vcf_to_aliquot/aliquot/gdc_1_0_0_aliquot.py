#!/usr/bin/env python3
"""gdc-1.0.0-aliquot"""

from maflib.column_types import MutationStatus

from aliquotmaf.vcf_to_aliquot.aliquot import Aliquot
from aliquotmaf.vcf_to_aliquot.annotators.cosmic import CosmicId
from aliquotmaf.vcf_to_aliquot.annotators.dbsnp_validation import DbSnpValidation
from aliquotmaf.vcf_to_aliquot.annotators.hotspot import Hotspot
from aliquotmaf.vcf_to_aliquot.annotators.nontcga_exac import NonTcgaExac
from aliquotmaf.vcf_to_aliquot.annotators.reference_context import ReferenceContext
from aliquotmaf.vcf_to_aliquot.filters.exac import ExAC
from aliquotmaf.vcf_to_aliquot.filters.gdc_blacklist import GdcBlacklist
from aliquotmaf.vcf_to_aliquot.filters.gdc_pon import GdcPon
from aliquotmaf.vcf_to_aliquot.filters.multiallelic import Multiallelic
from aliquotmaf.vcf_to_aliquot.filters.nonexonic import NonExonic
from aliquotmaf.vcf_to_aliquot.filters.normal_depth import NormalDepth
from aliquotmaf.vcf_to_aliquot.filters.offtarget import OffTarget


class GDC_1_0_0_Aliquot(Aliquot):

    VERSION = "gdc-1.0.0"
    ANNOTATION = "gdc-1.0.0-aliquot"

    ANNOTATORS = (
        MutationStatus,
        ReferenceContext,
        DbSnpValidation,
        CosmicId,
        NonTcgaExac,
        Hotspot,
    )
    FILTERS = (
        ExAC,
        GdcBlacklist,
        NormalDepth,
        GdcPon,
        Multiallelic,
        NonExonic,
        OffTarget,
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setup_annotators(self, args):
        """
        Sets up all annotator classes.
        """
        assert self.ANNOTATORS is not None, "Annotator classes must be defined"
        if not self.annotators:
            self.annotators = {}
            for Annotator in self.ANNOTATORS:
                self.annotators[Annotator.__name__] = Annotator.setup(self.scheme, args)

    def setup_filters(self, args):
        """
        Sets up all filter classes.
        """
        assert self.FILTERS is not None, "Filter classes must be defined"
        if not self.filters:
            self.filters = {}
            for Filter in self.FILTERS:
                filt_obj = Filter.setup(args)
                if filt_obj:
                    self.filters[Filter.__name__] = filt_obj


# __END__
