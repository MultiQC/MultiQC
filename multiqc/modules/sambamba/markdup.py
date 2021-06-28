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
        # Detect if sample is single or paired-end reads. Process differently.
        # Give user warning if redundant samples are found.

        self.markdup_data_dict = dict()

        for f in self.find_log_files("sambamba/markdup"):

            # define sample name
            sample_name = f["s_name"]

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
        Detect and calculate duplicate rate by single/paired-end samples.
        Outputs dictionary of markdup stats of a sample by single/paired-end type.
        This dict is associated with corresponding sample in key of markdup_data_dict.
        """

        """NOTE: reads sum to between 98 - 99% of input BAM file. Is Sambamba Markdup not recognizing the entire file?"""

        sorted_end_pairs = int(re.search("sorted \d+ end pairs", f).group(0).split(" ")[1])
        single_ends = int(re.search("and \d+ single ends", f).group(0).split(" ")[1])
        single_unmatched_pairs = int(re.search("among them \d+ unmatched", f).group(0).split(" ")[2])
        duplicate_reads = int(re.search("found \d+ duplicates", f).group(0).split(" ")[1])

        duplicate_rate = self.calculate_duplicate_rate(
            sortedEndPairs=sorted_end_pairs,
            singleEnds=single_ends,
            singleUnmatchedPairs=single_unmatched_pairs,
            duplicateReads=duplicate_reads,
        )

        sample_statistics = self.calculate_read_proportions(
            sortedEndPairs=sorted_end_pairs,
            singleEnds=single_ends,
            singleUnmatchedPairs=single_unmatched_pairs,
            duplicateReads=duplicate_reads,
            duplicateRate=duplicate_rate,
        )

        return sample_statistics

    def calculate_duplicate_rate(self, sortedEndPairs, singleEnds, singleUnmatchedPairs, duplicateReads):

        """
        Input all regexed stats from markdup log.
        Detect if sample is paired or single-end.
        Calculate duplicate rate as needed.
        Print debug msg if a sample neither fits singleEnds or PE scheme.
        Output duplicate rate.
        """

        if (sortedEndPairs > 0) and (singleUnmatchedPairs > 0):
            duplicate_rate = duplicateReads / (sortedEndPairs * 2 + singleEnds - singleUnmatchedPairs) * 100
            duplicate_rate = round(duplicate_rate, 2)
        elif (sortedEndPairs == 0) and (singleUnmatchedPairs == 0):
            duplicate_rate = duplicateReads / singleEnds * 100
            duplicate_rate = round(duplicate_rate, 2)
        else:
            log.debug("Sample {} is neither paired-end nor single-end.".format(f["s_name"]))
            duplicate_rate = 0

        return duplicate_rate

    def calculate_read_proportions(
        self, sortedEndPairs, singleEnds, singleUnmatchedPairs, duplicateReads, duplicateRate
    ):

        """
        Input all regexed stats from markdup log and duplicate rate.
        Detect if sample is PE or singleEnds.
        Return correct number of reads by type in context of showing proportions of a whole.
        """

        if (sortedEndPairs > 0) and (singleUnmatchedPairs > 0):
            stats = {
                "Total Sorted Paired End Reads": sortedEndPairs * 2 - duplicateReads,
                "Total Single End Reads": singleEnds - 2 * singleUnmatchedPairs,
                "Total Single Unmatched Reads": 2 * singleUnmatchedPairs,
                "Total Duplicate Reads": duplicateReads,
                "Duplicate Rates": duplicateRate,
            }
        elif (sortedEndPairs == 0) and (singleUnmatchedPairs == 0):
            stats = {
                "Total Sorted Paired End Reads": sortedEndPairs,  # 0
                "Total Single End Reads": singleEnds - duplicateReads,
                "Total Single Unmatched Reads": singleUnmatchedPairs,  # 0
                "Total Duplicate Reads": duplicateReads,
                "Duplicate Rates": duplicateRate,
            }
        else:
            stats = {
                "Total Sorted Paired End Reads": 0,
                "Total Single End Reads": 0,
                "Total Single Unmatched Reads": 0,
                "Total Duplicate Reads": 0,
                "Duplicate Rates": 0,
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
        headers["Total Duplicate Reads"] = {
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
        cats = [
            "Total Sorted Paired End Reads",
            "Total Single End Reads",
            "Total Single Unmatched Reads",
            "Total Duplicate Reads",
        ]

        config = {"id": "SambambaMarkdupBargraph", "title": "Sambamba Markdup: Duplicate Counts", "ylab": "# Reads"}

        self.add_section(
            name="Sambamba Markdup", anchor="SambambaMarkdup", plot=bargraph.plot(self.markdup_data_dict, cats, config)
        )
