#!/usr/bin/env python

"""MultiQC plugin to parse output from fadu"""

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import re

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="FADU",
            anchor="fadu",
            info="is a quantification tool designed specifically for prokaryotic RNA-Seq analyses.",
            href="https://igs.github.io/FADU/",
            doi="10.1128/mSystems.00917-20",
        )

        # Get all reports and load them
        self.count_data = dict()
        for f in self.find_log_files("fadu"):
            parsed_data = self.parse_fadu_data(f)
            if parsed_data is not None:
                s_name = f["s_name"]
                if s_name == "" or s_name == "Aligned.sortedByCoord.out.counts.txt":
                    s_name = self.clean_s_name(os.path.basename(f["root"]), f, root=os.path.dirname(f["root"]))
                if s_name in self.count_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.add_data_source(f)
            self.count_data[s_name] = parsed_data

        # Filter to strip out ignored sample names
        self.count_data = self.ignore_samples(self.count_data)

        if len(self.count_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.count_data)))

        if len(self.count_data) > 0:
            # Basic Stats Table
            self.fadu_stats_table()

            # Write parsed report data to a file
            self.write_data_file(self.count_data, "multiqc_fadu")

            # Plot number of detected genes
            self.add_section(
                name="Detected Genes",
                anchor="fadu_genes",
                description="number of detected genes with count > 0",
                plot=self.fadu_gene_chart(),
            )

    def parse_fadu_data(self, raw_data):
        """Parse the count file"""

        parsed_data = {}
        # data=raw_data

        # get number of aligments and counts for each gene:
        alig = list()
        counts = list()
        for line in raw_data["f"].splitlines():
            try:
                alig.append(float(line.split()[2]))
                counts.append(float(line.split()[3]))
            except:
                pass

        # Total number of BAM records that aligned to this feature's non-overlapping bases
        parsed_data["total_alig"] = sum(alig)

        # Number of detected genes (with count > 0)
        parsed_data["detected_genes"] = 0
        for val in counts:
            if val > 0:
                parsed_data["detected_genes"] += 1
        parsed_data["non_detected_genes"] = float(len(counts)) - parsed_data["detected_genes"]

        if len(parsed_data) == 0:
            return None
        return parsed_data

    def fadu_stats_table(self):
        """Take parsed data and add it to basic stats table"""

        headers = OrderedDict()
        headers["total_alig"] = {
            "title": "{} Assigned non-overlapping".format(config.read_count_prefix),
            "description": "Total Number of BAM records that aligned to features' non-overlapping bases ({})".format(
                config.read_count_desc
            ),
            "min": 0,
            "scale": "PuBu",
            "modify": lambda x: float(x) * config.read_count_multiplier,
            "shared_key": "read_count",
        }
        headers["detected_genes"] = {
            "title": "Detected genes",
            "description": "Detected genes with count > 0",
            "min": 0,
            "scale": "BuGn",
        }
        self.general_stats_addcols(self.count_data, headers)

    def fadu_gene_chart(self):
        """Make plot for detected genes"""

        # Specify data to depict
        keys = OrderedDict()
        keys["detected_genes"] = {"color": "#2988ff", "name": "Detected Genes"}
        keys["non_detected_genes"] = {"color": "#828181", "name": "Non-detected Genes"}

        # Config for plot
        pconfig = {
            "id": "fadu_detected_genes",
            "title": "FADU: Detected Genes (count > 0)",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Genes",
        }

        return bargraph.plot(self.count_data, keys, pconfig)
