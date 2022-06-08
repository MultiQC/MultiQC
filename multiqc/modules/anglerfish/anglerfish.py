#!/usr/bin/env python

""" MultiQC module to parse output from Anglerfish """

from __future__ import print_function
from collections import OrderedDict
import logging
import json

from multiqc import config
from multiqc.plots import bargraph, linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Anglerfish module class
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="anglerfish",
            anchor="anglerfish",
            href="https://github.com/remiolsen/anglerfish",
            info="A tool to assess Illumina libraries sequenced on Oxford Nanopore for the purpose of quality control.",
            doi=None,
        )

        # Find and load any anglerfish reports
        self.anglerfish_data = dict()
        self.anglerfish_pafstats = dict()
        self.anglerfish_samples = dict()
        self.anglerfish_undetermined = dict()

        for f in self.find_log_files("anglerfish", filehandles=True):
            self.parse_anglerfish_json(f)

        # Filter to strip out ignored sample names
        self.anglerfish_data = self.ignore_samples(self.anglerfish_data)

        # Stop execution of the data if no anglerfish data is found.
        if len(self.anglerfish_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.anglerfish_data)))

        # Write parsed report data to a file
        ## Parse whole JSON to save all its content
        self.write_data_file(self.anglerfish_data, "multiqc_anglerfish")

        # General Stats Table
        self.anglerfish_general_stats_table()

        # Paf Statistics bar plot
        self.add_section(
            name="Paf Statistics",
            anchor="anglerfish-paf-statistics",
            description="Paf Statistics of sampled reads.",
            plot=self.anglerfish_paf_stats_chart(),
        )
        # TODO: Sample statistics scatter plot
        self.add_section(
            name="Sample Statistics",
            anchor="anglerfish-sample-statistics",
            description="Outliers for the sample statistics",
            plot=self.anglerfish_sample_stats_chart(),
        )
        # TODO: Undetermined plot
        # self.add_section(
        #     name="Undetermined ??? TODO",
        #     anchor="TODO",
        #     description="TODO",
        #     plot=TODO,
        # )

    def parse_anglerfish_json(self, f):
        """Parse the JSON output from anglerfish and save the summary statistics"""
        try:
            parsed_json = json.load(f["f"])
        except:
            log.warning("Could not parse anglerfish JSON: '{}'".format(f["fn"]))
            return None

        # Fetch a sample name from the command
        s_name = f["s_name"]
        self.add_data_source(f, s_name)
        self.anglerfish_data[s_name] = {}
        self.anglerfish_pafstats[s_name] = {}
        self.anglerfish_samples[s_name] = {}
        self.anglerfish_undetermined[s_name] = {}
        ## TODO: Add function

    def anglerfish_general_stats_table(self):
        """TODO"""
        headers = OrderedDict()
        # TODO: Add headers

        self.general_stats_addcols(self.anglerfish_data, headers, "anglerfish")

    def anglerfish_paf_stats_chart(self):
        """TODO"""
        keys = OrderedDict()
        # TODO: Add keys
        # TODO: Add data structure for plot?
        # TODO: Add config
        return bargraph.plot(self.anglerfish_data, keys, config)

    def anglerfish_sample_stats_chart(self):
        """TODO"""
        keys = OrderedDict()
        # TODO: keys
        # TODO: ?? Maybe not, depending on plot: data structure for plot?
        # TODO: (p?)config
        # TODO: return type
