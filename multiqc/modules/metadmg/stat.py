# -*- coding: utf-8 -*-

""" MultiQC submodule to parse metadmg stat """

import logging
import pandas as pd

from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class StatReportMixin:
    """Mixin class, loaded by main metadmg MuliqcModule class."""

    def parse_metadmg_stat(self):
        """Find metadmg stat file and parse their data"""

        self.stat_types = ["nreads", "nalign", "mean_rlen", "var_rlen", "mean_gc", "var_gc"]

        self.metadmg_stat = dict()
        for f in self.find_log_files("metadmg/stat", filecontents=False, filehandles=False):
            # Read DF
            df = pd.read_table(
                f["fn"],
                sep="\t",
                usecols=["name", "rank"] + self.stat_types,
                index_col="name",
                comment=None,
            ).sort_values(by=self.sort_by, ascending=False)
            # Filter DF
            df = df[df["rank"] == self.rank].head(n=self.n_taxa).drop(columns="rank")

            parsed_data = df.to_dict()

            if len(parsed_data) > 0:
                if f["s_name"] in self.metadmg_stat:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))
                self.add_data_source(f, section="stat")
                self.metadmg_stat[f["s_name"]] = parsed_data

        # Filter to strip out ignored sample names
        self.metadmg_stat = self.ignore_samples(self.metadmg_stat)

        if len(self.metadmg_stat) == 0:
            return 0

        # Write parsed report data to a file
        self.write_data_file(self.metadmg_stat, "multiqc_metadmg_stat")

        # Add version
        self.add_software_version(None)

        # Make bargraph section
        self.barplot_section()

        # Return the number of logs that were found
        return len(self.metadmg_stat)

    def barplot_section(self):
        # Convert data
        stat_labels = {
            "nreads": "Number of reads",
            "nalign": "Number of alignments",
            "mean_rlen": "Mean read length",
            "var_rlen": "Read length variance",
            "mean_gc": "Mean GC content",
            "var_gc": "GC content variance",
        }

        data_plot = list()
        data_labels = list()
        for stat_type in self.stat_types:
            data_plot.append({s_name: data[stat_type] for s_name, data in self.metadmg_stat.items()})
            data_labels.append({"name": stat_labels[stat_type]})

        # Config for the plot
        pconfig = {
            "id": "metadmg_rank_plot",
            "hide_zero_cats": False,
            "title": "metaDMG: read statistics",
            "ylab": None,
            "use_legend": False,
            "tt_percentages": False,
            "data_labels": data_labels,
        }

        self.add_section(
            name="Read statistics by taxonomic rank",
            anchor="metadmg-rank",
            description=f"Read abundance statistics for top {self.n_taxa} {self.rank}",
            plot=bargraph.plot(data_plot, pconfig=pconfig),
        )
