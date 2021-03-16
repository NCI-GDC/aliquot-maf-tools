#!/usr/bin/env python3

from maflib.record import MafColumnRecord, MafRecord

from aliquotmaf.vcf_to_aliquot.aliquot.base import Aliquot
from aliquotmaf.vcf_to_aliquot.converters.builder import get_builder


def filter_maf(maf_record: MafRecord, aliquot: Aliquot) -> MafColumnRecord:
    gdc_filters = []
    for filt_obj in aliquot.filters.values():
        if filt_obj and filt_obj.filter(maf_record):
            gdc_filters.extend(filt_obj.tags)
    filter_record = get_builder(
        "GDC_FILTER", aliquot.scheme, value=";".join(sorted(gdc_filters))
    )
    return filter_record


# __END__
