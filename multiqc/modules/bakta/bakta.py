import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module analyses summary results from the Bakta annotation pipeline for bacterial genomes. The
    summary text file used is included in the Bakta output since v1.3.0. The MultiQC module was written for
    the output of v1.7.0.
    """

    def __init__(self):
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
            self.bakta[f["s_name"]] = self.parse_bakta(f)
            self.add_data_source(f)
        # Filter out ignored samples (given with --ignore-samples option)
        self.bakta = self.ignore_samples(self.bakta)

        # Raise error if dict is empty
        if len(self.bakta) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.bakta)} logs")

        # Write parsed report data to a file
        self.write_data_file(self.bakta, "multiqc_bakta")

        # Add most important Bakta annotation stats to general table
        headers = {
            "Count": {
                "title": "# contigs",
                "description": "Number of contigs",
                "min": 0,
                "scale": "Blues",
                "format": "{:,d}",
            },
            "Length": {
                "title": "# bases",
                "description": "Total number of bases in the contigs",
                "min": 0,
                "scale": "YlGn",
                "format": "{:,d}",
            },
            "CDSs": {
                "title": "# CDS",
                "description": "Number of found CDS",
                "min": 0,
                "scale": "YlGnBu",
                "format": "{:,d}",
                "shared_key": "gene_count",
            },
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
        for line in f["f"].splitlines():
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
