from aliquotmaf.subcommands.merge_aliquot.runners import GDC_1_0_0_Aliquot_Merged


class GDC_2_0_0_Aliquot_Merged(GDC_1_0_0_Aliquot_Merged):
    def __init__(self, options=dict()):
        super(GDC_2_0_0_Aliquot_Merged, self).__init__(options)

        # Schema
        self.options["version"] = "gdc-1.0.0"
        self.options["annotation"] = "gdc-2.0.0-aliquot-merged"

    @classmethod
    def __tool_name__(cls):
        return "gdc-2.0.0-aliquot-merged"
