# -*- coding: utf-8 -*-

""" MultiQC submodule to parse metadmg stat """

import logging
import pandas as pd

from multiqc import config
from multiqc.plots import bargraph, beeswarm

# Initialise the logger
log = logging.getLogger(__name__)



class StatReportMixin:
    """Mixin class, loaded by main metadmg MuliqcModule class."""

    def parse_metadmg_stat(self):
        """Find metadmg stat file and parse their data"""
   
        self.metadmg_stat = dict()
        for f in self.find_log_files("metadmg/stat", filecontents=False, filehandles=False):
            # Read DF
            df = pd.read_table(f["fn"], sep="\t", usecols=["name", "rank", "nalign", "nreads", "mean_rlen", "var_rlen", "mean_gc", "var_gc"], index_col=self.index, comment="#").sort_values(by=self.sort_by, ascending=False)
            # Filter DF
            df = df[df["rank"] == self.rank].head(n=self.n_taxa).drop(columns="rank")

            parsed_data = df.to_dict()
            log.debug(parsed_data)

            if len(parsed_data) > 0:
                if f["s_name"] in self.metadmg_stat:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))
                self.add_data_source(f, section="stat")
                self.metadmg_stat[f["s_name"]] = parsed_data
                
        # Filter to strip out ignored sample names
        self.metadmg_stat = self.ignore_samples(self.metadmg_stat)
        log.debug(self.metadmg_stat)

        if len(self.metadmg_stat) == 0:
            return 0

        # Write parsed report data to a file
        self.write_data_file(self.metadmg_stat, "multiqc_metadmg_stat")


        # Make dot-plot section
        self.dotplot_section()

        # Make bargraph section
        self.barplot_section()

        # Return the number of logs that were found
        return len(self.metadmg_stat)


    def dotplot_section(self):
        # Make dot plot of counts
        n_reads = {
            "min": 0,
            "suffix": "reads",
        }
        read_len = {
            "min": 0,
            "decimalPlaces": 4,
        }
        read_gc = {
            "min": 0,
            "max": 100,
            "decimalPlaces": 4,
        }
        keys = dict()
        keys["nreads"] = dict(n_reads, **{"title": "Total number of reads"})
        keys["nalign"] = dict(n_reads, **{"title": "Total number of alignments"})
        keys["mean_rlen"] = dict(read_len, **{"title": "Reads length mean"})
        keys["var_rlen"] = dict(read_len, **{"title": "Reads length variance"})
        keys["mean_gc"] = dict(read_gc, **{"title": "Reads GC mean"})
        keys["var_gc"] = dict(read_gc, **{"title": "Reads GC varaince"})

        self.add_section(
            name="Statistics",
            anchor="metadmg-stat",
            description="This module parses the output from <code>metadmg</code>.",
            plot=beeswarm.plot(self.metadmg_stat, keys, {"id": "metadmg-stat-dp"}),
        )


    def barplot_section(self):
        # Config for the plot
        pconfig = {
            "id": "metadmg_rank_plot",
            "title": f"metaDMG top {self.rank}",
            "ylab": self.rank.title(),
            "cpswitch_counts_label": "Number of Reads",
        }
                
        self.add_section(
            name="Rank distribution",
            anchor="metadmg-stat-rank",
            description=f"Top {self.n_taxa} {self.rank} by {self.sort_by}",
            plot=bargraph.plot(self.metadmg_stat[self.sort_by], pconfig=pconfig),
        )
