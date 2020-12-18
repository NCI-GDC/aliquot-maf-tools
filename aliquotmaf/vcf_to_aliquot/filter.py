#!/usr/bin/env python3

from aliquotmaf.vcf_to_aliquot.converters.builder import get_builder


def filter(maf_record: MafRecord, aliquot: Aliquot):
    gdc_filters = []
    for filt_key in aliquot.filters:
        filt_obj = aliquot.filters[filt_key]
        if filt_obj and filt_obj.filter(maf_record):
            gdc_filters.extend(filt_obj.tags)
    filter_record = get_builder(
        "GDC_FILTER", aliquot.scheme, value=";".join(sorted(gdc_filters))
    )


# __END__
