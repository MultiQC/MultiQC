#!/usr/bin/env python

from __future__ import print_function
from collections import OrderedDict
import os
import logging
import re

from multiqc import config
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class SambambaMarkdupMixin:

    """Find and parse Sambamba Markdup output log files and calculate duplication rate"""

    def parse_sambamba_markdup(self):

        # Clean sample name from 'markdup_sample_1' to 'sample_1'
        # Find and load sambamba logs to markdup_data_dict.
        # Regex for key phrases and calculate duplicate rate.
        # Give user warning if redundant samples are found.

        self.markdup_data_dict = dict()

        for f in self.find_log_files("sambamba/markdup"):

            # define sample name
            sample_name = self.clean_s_name(f["s_name"].replace("markdup_", ""), f["root"])

            # parse sambamba output by sample name
            self.markdup_data_dict[sample_name] = self.get_markdup_stats(f["f"])

            # filter away samples if MultiQC user does not want them
            self.markdup_data_dict = self.ignore_samples(self.markdup_data_dict)

            # add results to multiqc_sources.txt
            self.add_data_source(f)

            # warn user if duplicate samples found
            if f["s_name"] in self.markdup_data_dict:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))

        if len(self.markdup_data_dict) == 0:
            return 0

        # Write parsed files to a file
        self.markdup_write_data_to_file()

        # Add sambamba markdup to general statistics table
        self.markdup_general_stats_table()

        # Add sambamba markdup bargraph section
        self.markdup_section()

        # return the length of log files parsed
        return len(self.markdup_data_dict)

    def get_markdup_stats(self, f):

        """
        Input actual content of sambamba markdup log output.
        Regex for important phrases and extract reads.
        Calculate and return duplicate rate.
        Outputs dictionary of markdup stats of a sample.
        This dict is associated with corresponding sample in key of markdup_data_dict.
        """

        sorted_end_pairs = int(re.search("sorted \d+ end pairs", f).group(0).split(" ")[1])
        single_ends = int(re.search("and \d+ single ends", f).group(0).split(" ")[1])
        single_unmatched_pairs = int(re.search("among them \d+ unmatched", f).group(0).split(" ")[2])
        duplicate_reads = int(re.search("found \d+ duplicates", f).group(0).split(" ")[1])

        duplicate_rate = duplicate_reads / (sorted_end_pairs * 2 + single_ends - single_unmatched_pairs) * 100
        duplicate_rate = round(duplicate_rate, 2)

        stats = {
            "Sorted End Pairs": sorted_end_pairs,
            "Single Ends": single_ends,
            "Single Unmatched Pairs": single_unmatched_pairs,
            "Duplicate Reads": duplicate_reads,
            "Duplicate Rates": duplicate_rate,
        }

        return stats

    def markdup_write_data_to_file(self):

        """
        Write results to output file.
        """

        self.write_data_file(self.markdup_data_dict, "multiqc_markdup", sort_cols=True)

    def markdup_general_stats_table(self):

        """
        Take parsed stats from sambamba markdup to general stats table at the top of the report.
        """

        headers = OrderedDict()
        headers["Duplicate Reads"] = {
            "title": "Duplicate Reads",
            "description": "Number of Duplicate Reads per Sample",
            "scale": "RdYlGn-rev",
            "format": "{:,.0f}",
        }
        headers["Duplicate Rates"] = {
            "title": "Duplicate Rates",
            "description": "Rate of Duplication per Sample",
            "scale": "RdYlGn-rev",
            "format": "{:,.0f}",
            "suffix": "%",
        }
        self.general_stats_addcols(self.markdup_data_dict, headers)

    def markdup_section(self):

        """
        Add markdup statistics as bar graph to multiQC report.
        """

        # plot these categories, but not duplicate rate.
        cats = ["Sorted End Pairs", "Single Ends", "Single Unmatched Pairs", "Duplicate Reads"]

        config = {"id": "SambambaMarkdupBargraph", "title": "Sambamba Markdup: Duplicate Counts", "ylab": "# Reads"}

        self.add_section(
            name="Sambamba Markdup", anchor="SambambaMarkdup", plot=bargraph.plot(self.markdup_data_dict, cats, config)
        )
