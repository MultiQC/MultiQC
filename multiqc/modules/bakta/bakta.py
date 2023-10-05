""" MultiQC module to parse output from Bakta """


import logging
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
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
            info="is a tool for the rapid & standardized annotation of bacterial genomes, MAGs & plasmids",
            doi="10.1099/mgen.0.000685",
        )

        # Parse logs (txt files)
        self.bakta = dict()
        for f in self.find_log_files("bakta"):
            self.bakta[f["s_name"]] = self.parse_bakta(f)
            self.add_data_source(f)
        # Filter out ignored samples (given with --ignore-samples option)
        self.bakta = self.ignore_samples(self.bakta)

        # Raise error if dict is empty
        if len(self.bakta) == 0:
            raise ModuleNoSamplesFound

        log.info("Found {} logs".format(len(self.bakta)))

        # Write parsed report data to a file
        self.write_data_file(self.bakta, "multiqc_bakta")

        # Add most important Bakta annotation stats to general table
        headers = OrderedDict()
        headers["Count"] = {
            "title": "# contigs",
            "description": "Number of contigs",
            "min": 0,
            "scale": "Blues",
            "format": "{:,d}",
        }
        headers["Length"] = {
            "title": "# bases",
            "description": "Total number of bases in the contigs",
            "min": 0,
            "scale": "YlGn",
            "format": "{:,d}",
        }
        headers["CDSs"] = {
            "title": "# CDS",
            "description": "Number of found CDS",
            "min": 0,
            "scale": "YlGnBu",
            "format": "{:,d}",
            "shared_key": "gene_count",
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

        metrics = (
            "Length",
            "Count",
            "N50",
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

        data = {}
        for line in f["contents_lines"]:
            s = line.strip().split(": ")
            if s[0] in metrics:
                data[s[0].replace(" ", "_")] = int(s[1])

            if s[0] == "Software":
                self.add_software_version(s[1].strip("v"), f["s_name"])
        return data

    def bakta_barplot(self):
        """Make a basic plot of the annotation stats"""

        # Specify the order of the different categories
        categories = {
            metric: {"name": cat}
            for metric, cat in {
                "CDSs": "CDS",
                "rRNAs": "rRNA",
                "tRNAs": "tRNA",
                "tmRNAs": "tmRNA",
                "ncRNAs": "ncRNA",
                "ncRNA_regions": "ncRNA region",
                "hypotheticals": "hypothetical",
                "CRISPR_arrays": "CRISPR array",
                "pseudogenes": "pseudogene",
                "signal_peptides": "signal peptide",
                "sORFs": "sORF",
                "gaps": "gap",
                "oriCs": "oriC",
                "oriVs": "oriV",
                "oriTs": "oriT",
            }.items()
        }

        plot_config = {
            "id": "bakta_feature_types",
            "title": "Bakta: Feature Types",
            "ylab": "# Contigs",
            "cpswitch_counts_label": "Features",
        }
        return bargraph.plot(self.bakta, categories, plot_config)
