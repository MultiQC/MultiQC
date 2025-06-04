from multiqc.base_module import BaseMultiqcModule

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="BD csv",
            anchor="bd_custom_csv",
            info="Parses a CSV file to extract summary metrics and display them in MultiQC.",
        )