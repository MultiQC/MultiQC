import logging
import os
import re

import numpy as np

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="JCVI Genome Annotation",
            anchor="jcvi",
            href="https://pypi.org/project/jcvi/",
            info="computes statistics on genome annotation.",
            doi="10.5281/zenodo.31631",
        )

        # Parse logs
        self.jcvi = dict()
        for f in self.find_log_files("jcvi", filehandles=True):
            self.parse_jcvi(f)

        # Filter to strip out ignored sample names
        self.jcvi = self.ignore_samples(self.jcvi)

        if len(self.jcvi) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.jcvi)} logs")
        self.write_data_file(self.jcvi, "multiqc_jcvi")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Add most important JCVI annotation stats to the general table
        headers = {
            "genes": {
                "title": "Number of genes",
                "description": "Number of genes",
                "format": "{:i}",
            },
            "transcripts": {
                "title": "Number of transcripts",
                "description": "Number of transcripts",
                "format": "{:i}",
            },
            "mean_gene_size": {
                "title": "Mean gene length (bp)",
                "description": "Mean gene length",
                "format": "{:i}",
            },
        }
        self.general_stats_addcols(self.jcvi, headers)

        self.add_section(
            plot=self.jcvi_barplot_feature_counts(),
            name="Number of features",
            description="Total number of genes, transcripts and exons found in each dataset.",
        )

        self.add_section(
            plot=self.jcvi_barplot_feature_lengths(),
            name="Mean size of features",
            description="Mean size of all genes, transcripts and exons found in each dataset.",
        )

        self.add_section(
            plot=self.jcvi_barplot_features_per_gene(),
            name="Features per gene",
            description="Mean and maximum number of transcripts and exons per gene.",
        )

        self.add_section(
            plot=self.jcvi_barplot_isoforms(),
            name="Isoforms",
            description="Number of genes found to have multiple isoforms",
        )

        feature_length_plot = self.jcvi_linegraph_feature_length()
        if feature_length_plot:
            self.add_section(
                name="Feature length distribution",
                anchor="jcvi_feature_length",
                description="This plot shows the distribution of gene, exon and intron length (if available).",
                helptext="""
                Values are binned in buckets as follows:

                * Gene length: 100bp (e.g. the value at 150bp represents the number of genes having a length between 100 and 199bp).
                * Exon length: 25bp (e.g. the value at 112bp represents the number of exons having a length between 100 and 124bp).
                * Intron length: 25bp (e.g. the value at 112bp represents the number of introns having a length between 100 and 124bp).
                """,
                plot=feature_length_plot,
            )

        exon_count_plot = self.jcvi_linegraph_exon_count()
        if exon_count_plot:
            self.add_section(
                name="Exon count distribution",
                anchor="jcvi_exon_count",
                description="This plot shows the distribution of exon number per genes.",
                helptext="""
                If you look at y-axis value corresponding to the x-axis value "3", you get the number of genes having 3 exons.
                """,
                plot=exon_count_plot,
            )

    def parse_jcvi(self, f):
        # Look at the first three lines, they are always the same
        first_line = f["f"].readline()
        second_line = f["f"].readline()
        third_line = f["f"].readline()

        # If any of these fail, it's probably not a jcvi summary file
        if not all(
            (
                first_line.startswith("==================================================================="),
                second_line.startswith("                                                   o            all"),
                third_line.startswith("-------------------------------------------------------------------"),
            )
        ):
            return

        # Set up sample dict
        s_name = f["s_name"]
        if s_name in self.jcvi:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        self.jcvi[s_name] = dict()

        # Define parsing regexes
        regexes = {
            "transcripts_per_gene": r"Max number of transcripts per gene\s+([\d,]+)",
            "mean_exons_size": r"Mean exon size\s+([\d,.]+)",
            "mean_gene_size": r"Mean gene locus size \(first to last exon\)\s+([\d,.]+)",
            "mean_exon_number": r"Mean number of distinct exons per gene\s+([\d,.]+)",
            "mean_transcript_number": r"Mean number of transcripts per gene\s+([\d,.]+)",
            "mean_transcript_size": r"Mean transcript size \(UTR, CDS\)\s+([\d,.]+)",
            "exons": r"Number of distinct exons\s+([\d,]+)",
            "genes": r"Number of genes\s+([\d,]+)",
            "genes_with_alt": r"Number of genes with alternative transcript variants\s+(\d+)",
            "genes_with_alt_percent": r"Number of genes with alternative transcript variants\s+\d+\s\(([\d.]+)%\)",
            "multiexon_genes": r"Number of multi-exon genes\s+([\d,]+)",
            "multiexon_genes_percent": r"Number of multi-exon genes\s+[\d,]+\s\(([\d.]+)%\)",
            "transcripts": r"Number of predicted transcripts\s+([\d,]+)",
            "singleexon_genes": r"Number of single-exon genes\s+([\d,]+)",
            "singleexon_genes_percent": r"Number of single-exon genes\s+[\d,]+\s\(([\d.]+)%\)",
        }

        # Go through lines in the file checking regexes
        for line in f["f"]:
            for key in list(regexes.keys()):
                match = re.search(regexes[key], line)
                if match:
                    try:
                        self.jcvi[s_name][key] = int(match.group(1).replace(",", ""))
                    except ValueError:
                        self.jcvi[s_name][key] = float(match.group(1).replace(",", ""))

                    # Don't keep searching for this regex in this file, now that we've found it
                    del regexes[key]
                    break

        # Parse histograms
        hist_configs = {
            "exon_count": ("Exon_Count", 1),
            "exon_length": ("Exon_Length", 25),
            "gene_length": ("Gene_Length", 100),
            "intron_length": ("Intron_Length", 25),
        }
        for key, config in hist_configs.items():
            results = self.parse_hists(os.path.join(f["root"], config[0], s_name + ".txt"), bin_by=config[1])
            if results:
                self.jcvi[s_name][key] = results

        # Check that we managed to parse some data for this sample
        if len(self.jcvi[s_name]) == 0:
            del self.jcvi[s_name]
            return

        self.add_data_source(f, s_name)

    def parse_hists(self, stat_file, bin_by=1):
        stat_table = {}

        vals = []
        if os.path.isfile(stat_file):
            with open(stat_file, "r") as stat_handle:
                for stat_line in stat_handle:
                    vals.append(int(stat_line))

        if not vals:
            return stat_table

        # Generate bin list
        bins = list(range(1, max(vals) + bin_by, bin_by))

        # Assign each value to its bin
        inds = np.digitize(vals, bins, right=True)

        for ind in inds:
            val = bins[int(ind)]
            if val not in stat_table:
                stat_table[val] = 1
            else:
                stat_table[val] += 1

        # Add missing values
        if stat_table:
            for required_val in bins:
                if required_val not in stat_table:
                    stat_table[required_val] = 0

        return stat_table

    def jcvi_barplot_feature_counts(self):
        plot_config = {
            "id": "jcvi_plot_feature_counts_plot",
            "title": "JCVI: Number of features",
            "ylab": "Number of Genes",
            "data_labels": [
                {"name": "Genes", "ylab": "Number of Genes"},
                {"name": "Transcripts", "ylab": "Number of Transcripts"},
                {"name": "Exons", "ylab": "Number of Exons"},
            ],
        }
        cats = [
            {
                "multiexon_genes": {"name": "Multi-exon"},
                "singleexon_genes": {"name": "Single-exon"},
            },
            {"transcripts": {"name": "Transcripts"}},
            {"exons": {"name": "Exons"}},
        ]
        return bargraph.plot([self.jcvi, self.jcvi, self.jcvi], cats, plot_config)

    def jcvi_barplot_feature_lengths(self):
        plot_config = {
            "id": "jcvi_plot_features_len",
            "title": "JCVI: Mean sizes of features",
            "ylab": "Base pairs",
            "cpswitch": False,
            "tt_decimals": 1,
            "data_labels": [
                {"name": "Genes", "ylab": "Base pairs"},
                {"name": "Transcripts", "ylab": "Base pairs"},
                {"name": "Exons", "ylab": "Base pairs"},
            ],
        }
        cats = [
            {"mean_gene_size": {"name": "Mean gene size"}},
            {"mean_transcript_size": {"name": "Mean transcript size"}},
            {"mean_exons_size": {"name": "Mean exons size"}},
        ]
        return bargraph.plot([self.jcvi, self.jcvi, self.jcvi], cats, plot_config)

    def jcvi_barplot_features_per_gene(self):
        cats = [
            {
                "mean_transcript_number": {"name": "Mean transcripts per gene"},
                "transcripts_per_gene": {"name": "Maximum transcripts per gene"},
            },
            {"mean_exon_number": {"name": "Mean exons per genes"}},
        ]
        plot_config = {
            "id": "jcvi_plot_features_per_genes",
            "title": "JCVI: Features per gene",
            "ylab": "# Transcripts per gene",
            "cpswitch": False,
            "tt_decimals": 1,
            "data_labels": [
                {"name": "Transcripts", "ylab": "# Transcripts per gene"},
                {"name": "Exons", "ylab": "# Exons per gene"},
            ],
        }

        return bargraph.plot([self.jcvi, self.jcvi], cats, plot_config)

    def jcvi_barplot_isoforms(self):
        keys = {"genes_with_alt": {"name": "Genes with multiple isoforms"}}

        plot_config = {
            "id": "jcvi_plot_isoforms",
            "title": "JCVI: Genes with multiple isoforms",
            "ylab": "# genes with multiple isoforms",
            "cpswitch": False,
        }

        return bargraph.plot(self.jcvi, keys, plot_config)

    def jcvi_linegraph_feature_length(self):
        plot_config = {
            "id": "jcvi_feature_length_plot",
            "title": "JCVI: Feature length repartition",
            "ylab": "# Genes",
            "xlab": "Gene length",
            "xsuffix": "bp",
            "data_labels": [
                {"name": "Genes", "ylab": "# Genes", "xlab": "Gene length"},
                {"name": "Exons", "ylab": "# Exons", "xlab": "Exon length"},
                {"name": "Introns", "ylab": "# Introns", "xlab": "Intron length"},
            ],
        }

        genes_plot_data = {x: self.jcvi[x]["gene_length"] for x in self.jcvi if "gene_length" in self.jcvi[x]}
        exon_plot_data = {x: self.jcvi[x]["exon_length"] for x in self.jcvi if "exon_length" in self.jcvi[x]}
        intron_plot_data = {x: self.jcvi[x]["intron_length"] for x in self.jcvi if "intron_length" in self.jcvi[x]}

        plot_data = []
        if genes_plot_data:
            plot_data.append(genes_plot_data)
        if exon_plot_data:
            plot_data.append(exon_plot_data)
        if intron_plot_data:
            plot_data.append(intron_plot_data)

        if len(plot_data) == 0:
            return None

        return linegraph.plot(plot_data, plot_config)

    def jcvi_linegraph_exon_count(self):
        plot_config = {
            "id": "jcvi_exon_count_plot",
            "title": "JCVI: Exon count repartition",
            "ylab": "# genes",
            "xlab": "Exon count",
            "xDecimals": False,
            "ymin": 0,
        }

        plot_data = {x: self.jcvi[x]["exon_count"] for x in self.jcvi if "exon_count" in self.jcvi[x]}

        if not plot_data:
            return None

        return linegraph.plot(plot_data, plot_config)
