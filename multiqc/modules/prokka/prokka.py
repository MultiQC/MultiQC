""" MultiQC module to parse output from Prokka """


import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Prokka",
            anchor="prokka",
            href="http://www.vicbioinformatics.com/software.prokka.shtml",
            info="is a software tool for the rapid annotation of prokaryotic genomes.",
            doi="10.1093/bioinformatics/btu153",
        )

        # Parse logs
        self.prokka = dict()
        for f in self.find_log_files("prokka", filehandles=True):
            self.parse_prokka(f)

        # Filter to strip out ignored sample names
        self.prokka = self.ignore_samples(self.prokka)

        if len(self.prokka) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.prokka)} logs")

        self.write_data_file(self.prokka, "multiqc_prokka")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Add most important Prokka annotation stats to the general table
        headers = {
            "organism": {
                "title": "Organism",
                "description": "Organism",
            },
            "contigs": {
                "title": "Contigs",
                "description": "Number of contigs",
                "min": 0,
            },
            "bases": {
                "title": "Bases",
                "description": "Number of bases",
                "min": 0,
                "format": "{:i}%",
                "hidden": True,
            },
            "CDS": {
                "title": "CDS",
                "description": "Number of CDS",
                "min": 0,
                "format": "{:i}%",
            },
        }
        self.general_stats_addcols(self.prokka, headers)

        # User can set configuration attributes, 'prokka_table', and
        # 'prokka_barplot', to specify whether to include a table or a barplot, or both.
        # Default is to make a plot only.
        if getattr(config, "prokka_table", False):
            self.add_section(plot=self.prokka_table())
        if getattr(config, "prokka_barplot", True):
            descr_plot = "This barplot shows the distribution of different types of features found in each contig."
            helptext = """
            `Prokka` can detect different features:

            - CDS
            - rRNA
            - tmRNA
            - tRNA
            - miscRNA
            - signal peptides
            - CRISPR arrays

            This barplot shows you the distribution of these different types of features found in each contig.
            """
            self.add_section(plot=self.prokka_barplot(), helptext=helptext, description=descr_plot)

    def parse_prokka(self, f):
        """Parse prokka txt summary files.

        Prokka summary files are difficult to identify as there are practically
        no distinct prokka identifiers in the filenames or file contents. This
        parser makes an attempt using the first three lines, expected to contain
        organism, contigs, and bases statistics.
        """

        s_name = None

        # Look at the first three lines, they are always the same
        first_line = f["f"].readline()
        contigs_line = f["f"].readline()
        bases_line = f["f"].readline()
        # If any of these fail, it's probably not a prokka summary file
        if not all(
            (first_line.startswith("organism:"), contigs_line.startswith("contigs:"), bases_line.startswith("bases:"))
        ):
            return

        # Get organism and sample name from the first line
        # Assumes organism name only consists of two words,
        # i.e. 'Genusname speciesname', and that the remaining
        # text on the organism line is the sample name.
        try:
            organism = " ".join(first_line.strip().split(":", 1)[1].split()[:2])
            s_name = self.clean_s_name(" ".join(first_line.split()[3:]), f)
        except KeyError:
            organism = first_line.strip().split(":", 1)[1]
            s_name = f["s_name"]
        # Don't try to guess sample name if requested in the config
        if getattr(config, "prokka_fn_snames", False):
            s_name = f["s_name"]

        if s_name in self.prokka:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        self.prokka[s_name] = dict()
        self.prokka[s_name]["organism"] = organism
        self.prokka[s_name]["contigs"] = int(contigs_line.split(":")[1])
        self.prokka[s_name]["bases"] = int(bases_line.split(":")[1])

        # Get additional info from remaining lines
        for line in f["f"]:
            description, value = line.split(":")
            try:
                self.prokka[s_name][description] = int(value)
            except ValueError:
                log.warning("Unable to parse line: '%s'", line)

        self.add_data_source(f, s_name)

    def prokka_table(self):
        """Make basic table of the annotation stats"""

        # Specify the order of the different possible categories
        headers = {
            "organism": {
                "title": "Organism",
                "description": "Organism name",
            },
            "contigs": {
                "title": "# contigs",
                "description": "Number of contigs in assembly",
                "format": "{:i}",
            },
            "bases": {
                "title": "# bases",
                "description": "Number of nucleotide bases in assembly",
                "format": "{:i}",
            },
            "CDS": {
                "title": "# CDS",
                "description": "Number of annotated CDS",
                "format": "{:i}",
            },
            "rRNA": {
                "title": "# rRNA",
                "description": "Number of annotated rRNA",
                "format": "{:i}",
            },
            "tRNA": {
                "title": "# tRNA",
                "description": "Number of annotated tRNA",
                "format": "{:i}",
            },
            "tmRNA": {
                "title": "# tmRNA",
                "description": "Number of annotated tmRNA",
                "format": "{:i}",
            },
            "misc_RNA": {
                "title": "# misc RNA",
                "description": "Number of annotated misc. RNA",
                "format": "{:i}",
            },
            "sig_peptide": {
                "title": "# sig_peptide",
                "description": "Number of annotated sig_peptide",
                "format": "{:i}",
            },
            "repeat_region": {
                "title": "# CRISPR arrays",
                "description": "Number of annotated CRSIPR arrays",
                "format": "{:i}",
            },
        }
        table_config = {
            "namespace": "prokka",
            "min": 0,
        }

        return table.plot(self.prokka, headers, table_config)

    def prokka_barplot(self):
        """Make a basic plot of the annotation stats"""

        # Specify the order of the different categories
        keys = {
            "CDS": {"name": "CDS"},
            "rRNA": {"name": "rRNA"},
            "tRNA": {"name": "tRNA"},
            "tmRNA": {"name": "tmRNA"},
            "misc_RNA": {"name": "misc RNA"},
            "sig_peptide": {"name": "Signal peptides"},
            "repeat_region": {"name": "CRISPR array"},
        }

        plot_config = {
            "id": "prokka_plot",
            "title": "Prokka: Feature Types",
            "ylab": "# Counts",
            "cpswitch_counts_label": "Features",
        }

        return bargraph.plot(self.prokka, keys, plot_config)
