#!/usr/bin/env python3
"""Main vcf2maf logic for spec gdc-1.0.0-aliquot"""

import pysam
from maflib.column_types import MutationStatus

from aliquotmaf.utils.utils import (
    assert_sample_in_header,
    extract_annotation_from_header,
)
from aliquotmaf.vcf_to_aliquot.aliquot import Aliquot
from aliquotmaf.vcf_to_aliquot.annotators.cosmic import CosmicId
from aliquotmaf.vcf_to_aliquot.annotators.dbsnp_validation import DbSnpValidation
from aliquotmaf.vcf_to_aliquot.annotators.hotspot import Hotspot
from aliquotmaf.vcf_to_aliquot.annotators.nontcga_exac import NonTcgaExac
from aliquotmaf.vcf_to_aliquot.annotators.reference_context import ReferenceContext
from aliquotmaf.vcf_to_aliquot.converters.utils import get_columns_from_header
from aliquotmaf.vcf_to_aliquot.filters import (
    ExAC,
    GdcBlacklist,
    GdcPon,
    Multiallelic,
    NonExonic,
    NormalDepth,
    OffTarget,
)


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
        if not self.annotators:
            self.annotators = {}
            for Annotator in self.ANNOTATORS:
                self.annotators[Annotator.__name__] = Annotator.setup(self.scheme, args)

    def setup_filters(self, args):
        """
        Sets up all filter classes.
        """
        if not self.filters:
            self.filters = {}
            for Filter in self.FILTERS:
                filt = Filter.setup(args)
                if filt:
                    self.filters[Filter.__name__] = filt


# __END__
