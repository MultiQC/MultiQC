""" MultiQC module to parse output from Bakta """


import logging
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Bakta",
            anchor="bakta",
            href="https://github.com/oschwengers/bakta",
            info="Rapid & standardized annotation of bacterial genomes, MAGs & plasmids",
            doi="10.1099/mgen.0.000685",
        )

        # Parse logs (txt files)
        self.bakta = dict()
        for f in self.find_log_files("bakta"):
            self.bakta[f["s_name"]] = self.parse_bakta(f["f"])
            self.add_data_source(f)
        # Filter out ignored samples (given with --ignore-samples option)
        self.bakta = self.ignore_samples(self.bakta)

        # Raise error if dict is empty
        if len(self.bakta) == 0:
            raise UserWarning

        log.info("Found {} logs".format(len(self.bakta)))

        # Write parsed report data to a file
        self.write_data_file(self.bakta, "multiqc_bakta")

        # Add most important Bakta annotation stats to general table
        headers = OrderedDict()
        headers["Count"] = {
            "title": "# contigs",
            "description": "Number of contigs",
            "min": 0,
            "format": "{:i}%",
        }
        headers["Length"] = {
            "title": "# bases",
            "description": "Number of bases",
            "min": 0,
            "format": "{:i}%",
        }
        headers["CDSs"] = {
            "title": "# CDS",
            "description": "Number of CDS",
            "min": 0,
            "format": "{:i}%",
        }
        self.general_stats_addcols(self.bakta, headers)

        # Helptext for Bakta barplot:
        descr_plot = "This barplot shows the distribution of different types of features found in each sample."

        helptext = """
           `Bakta` can detect different features:

            - tRNAs
            - tmRNAs
            - rRNAs
            - ncRNAs
            - ncRNA regions
            - CRISPR arrays
            - CDSs
            - pseudogenes
            - hypotheticals
            - signal peptides
            - sORFs
            - gaps
            - oriCs
            - oriVs
            - oriTs


            This barplot shows you the distribution of these different types of features found in each sample.
            """
        self.add_section(plot=self.bakta_barplot(), helptext=helptext, description=descr_plot)

    def parse_bakta(self, f):
        """Parse bakta txt summary files for annotation parameters."""

        params = (
            "Length",
            "Count",
            "tRNAs",
            "tmRNAs",
            "rRNAs",
            "ncRNAs",
            "ncRNA regions",
            "CRISPR arrays",
            "CDSs",
            "pseudogenes",
            "hypotheticals",
            "signal peptides",
            "sORFs",
            "gaps",
            "oriCs",
            "oriVs",
            "oriTs",
        )

        bakta_params = {}

        for l in f.splitlines():
            s = l.split(": ")
            if s[0] in params:
                bakta_params[s[0]] = int(s[1])
        return bakta_params

    def bakta_barplot(self):
        """Make a basic plot of the annotation stats"""

        # Specify the order of the different categories
        keys = OrderedDict()
        keys["CDSs"] = {"name": "CDS"}
        keys["rRNAs"] = {"name": "rRNA"}
        keys["tRNAs"] = {"name": "tRNA"}
        keys["tmRNAs"] = {"name": "tmRNA"}
        keys["ncRNAs"] = {"name": "ncRNA"}
        keys["ncRNA regions"] = {"name": "ncRNA region"}
        keys["hypotheticals"] = {"name": "hypothetical"}
        keys["CRISPR arrays"] = {"name": "CRISPR array"}
        keys["pseudogenes"] = {"name": "pseudogene"}
        keys["signal peptides"] = {"name": "signal peptide"}
        keys["sORFs"] = {"name": "sORF"}
        keys["gaps"] = {"name": "gap"}
        keys["oriCs"] = {"name": "oriC"}
        keys["oriVs"] = {"name": "oriV"}
        keys["oriTs"] = {"name": "oriT"}

        plot_config = {
            "id": "bakta_plot",
            "title": "Bakta: Feature Types",
            "ylab": " Counts",
            "cpswitch_counts_label": "Features",
        }
        return bargraph.plot(self.bakta, keys, plot_config)
