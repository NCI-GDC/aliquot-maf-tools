"""
Implements the dbSNP validation status annotations.
"""
from __future__ import absolute_import

import sqlite3 as lite

from .annotator import Annotator

from aliquotmaf.converters.builder import get_builder


class DbSnpValidation(Annotator):
    def __init__(self, scheme, source):
        super().__init__(name="DbSnpValidation", source=source, scheme=scheme)
        self.conn = None
        self.cur = None

    @classmethod
    def setup(cls, scheme, source):
        curr = cls(scheme, source)
        curr.logger.info("Connecting to dbsnp priority DB")
        curr.conn = lite.connect(source)
        curr.cur = curr.conn.cursor()
        return curr

    def annotate(self, maf_record):
        validation = None
        dbsnp_rs = maf_record["dbSNP_RS"].value
        if dbsnp_rs and dbsnp_rs != ["novel"]:
            results = []
            for dbsnp in dbsnp_rs:
                self.cur.execute(
                    "SELECT valstatus FROM dbsnpvalstat WHERE rsid=?", (dbsnp,)
                )
                for row in self.cur.fetchall():
                    results.append(row)
            if results:
                validation = ";".join([i[0] for i in sorted(list(set(results)))])
        maf_record["dbSNP_Val_Status"] = get_builder(
            "dbSNP_Val_Status", self.scheme, value=validation
        )
        return maf_record

    def shutdown(self):
        self.conn.close()
