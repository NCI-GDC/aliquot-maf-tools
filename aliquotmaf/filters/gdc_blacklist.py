"""
Applies the GDC Blacklist filter.
"""

from __future__ import absolute_import

import gzip

from .filter_base import Filter


class GdcBlacklist(Filter):
    def __init__(self, source, data):
        super().__init__(name="GDCBlacklist", source=source)
        self.tags = []
        self.data = data
        self.logger.info("Using GDC Blacklist {0}".format(source))

    @classmethod
    def setup(cls, source):
        # Load blacklist
        data = {}
        head = []
        reader = gzip.open if source.endswith(".gz") else open
        linect = 0
        for line in reader(source, "rt"):
            linect += 1
            if not head:
                head = line.rstrip().lower().split("\t")
                assert "tumor_aliquot_id" in head, (
                    'Required column "tumor_aliquot_id" missing from blacklist file {0}'.format(
                        cls.source
                    )
                )
                assert "tag" in head, (
                    'Required column "tag" missing from blacklist file {0}'.format(
                        cls.source
                    )
                )
            else:
                # Parse row into dict
                row = dict(zip(head, line.rstrip("\r\n").split("\t")))
                aliquot = row["tumor_aliquot_id"]
                tags = row["tag"].split(";")
                assert aliquot, (
                    "Missing required field tumor_aliquot_id in blacklist file - Line {0}".format(
                        linect
                    )
                )
                # Set dict tumor_aliquot_id -> tag
                if tags:
                    data[aliquot] = tags

        curr = cls(source, data)
        return curr

    def filter(self, maf_record):
        self.tags = []
        flag = False
        tumor_aliquot = str(maf_record["Tumor_Sample_UUID"].value)
        if tumor_aliquot in self.data:
            self.tags = self.data[tumor_aliquot]
            flag = True
        return flag

    def shutdown(self):
        pass
