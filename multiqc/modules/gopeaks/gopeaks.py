#!/usr/bin/env python

""" MultiQC module to parse output from gopeaks """

from __future__ import print_function
from collections import OrderedDict
import os
import logging
import json

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """gopeaks module"""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="GoPeaks",
            anchor="gopeaks",
            href="https://github.com/maxsonBraunLab/gopeaks",
            info="is designed to call peaks from aligned CUT&TAG sequencing reads.\
                  Gopeaks uses a binomial distribution to model the read counts \
                  in sliding windows across the genome and calculate peak regions \
                  that are enriched over the background.",
            doi=["10.1101/2022.01.10.475735"],
        )

        # data vars ---------------------------------------------------------------------

        # Find and load any gopeaks reports
        self.gopeaks_data = dict()
        for f in self.find_log_files("gopeaks", filehandles=True):

            parsed = self.parse_gopeaks_log(f)

            if parsed is not None:
                self.gopeaks_data[f["s_name"]] = parsed

            # filter away samples if MultiQC user does not want them
            self.gopeaks_data = self.ignore_samples(self.gopeaks_data)

            if f["s_name"] in self.gopeaks_data:
                log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f["fn"], f["s_name"]))

        if len(self.gopeaks_data) == 0:
            raise UserWarning
        if len(self.gopeaks_data) > 0:
            log.info("Found {} samples".format(len(self.gopeaks_data)))

        # Add sample log info to basic stats table
        self.gopeaks_general_stats_table()

        # Add sample log info to bargraph
        self.gopeaks_bargraph()

    # parsing functions -------------------------------------------------------------

    def parse_gopeaks_log(self, f):

        """
        Read gopeaks json log file and extract number of peaks.
        """

        d = {}
        path_to_file = f["root"] + "/" + f["fn"]
        f = open(path_to_file)
        sample_log = json.load(f)
        d["peak_counts"] = sample_log["peak_counts"]

        return d

    def gopeaks_general_stats_table(self):

        """
        Put peak counts to the general table.
        """

        headers = OrderedDict()
        headers["peak_counts"] = {
            "title": "GoPeaks Peak Counts",
            "description": "Number of peaks per sample",
            "min": 0,
            "scale": "YlGnBu",
            "format": "{:,.0f}",
        }
        self.general_stats_addcols(self.gopeaks_data, headers)

    def gopeaks_bargraph(self):

        """
        Put peak counts to a bargraph.
        """

        cats = OrderedDict()
        cats["peak_counts"] = {"name": "Peak Counts"}

        config = {
            "id": "GoPeaksBarGraph",
            "title": "GoPeaks: Number of Peaks by Sample",
            "ylab": "Sample",
            "logswitch": True,
            "cpswitch": False,
        }

        self.add_section(
            name="GoPeaks",
            anchor="gopeaks_bargraph",
            description="Number of peaks called by GoPeaks",
            plot=bargraph.plot(self.gopeaks_data, cats, config),
        )
