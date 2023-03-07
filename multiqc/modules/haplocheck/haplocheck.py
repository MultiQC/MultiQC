from multiqc.modules.base_module import BaseMultiqcModule


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Haplocheck",
            anchor="haplocheck",
            href="https://github.com/genepi/haplocheck",
            info="Haplocheck detects contamination patterns in mtDNA AND WGS sequencing studies by analyzing the mitochondrial DNA.",
            doi="10.1101/gr.256545.119",
        )
