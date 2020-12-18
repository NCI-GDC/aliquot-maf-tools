"""
Implements the dbSNP validation status annotations.
"""

import sqlite3

from aliquotmaf.annotators.annotator import Annotator
from aliquotmaf.converters.builder import get_builder


class DbSnpValidation(Annotator):
    def __init__(self, scheme, source):
        super().__init__(name="DbSnpValidation", source=source, scheme=scheme)
        self.conn = sqlite3.connect(source)
        self.cur = self.conn.cursor()

    @classmethod
    def setup(cls, scheme, args):
        curr = cls(scheme, args.dbsnp_priority_db)
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
        dbsnp_val_status_record = get_builder(
            "dbSNP_Val_Status", self.scheme, value=validation
        )
        return dbsnp_val_status_record

    def shutdown(self):
        self.conn.close()