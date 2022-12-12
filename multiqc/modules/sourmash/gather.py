#!/usr/bin/env python

""" MultiQC module to parse similarity matrix output by sourmash gather """

import csv
import logging
from collections import OrderedDict

from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class gather:
    def parse_gather(self):
        """
        Modeled after the kraken module and the sourmash compare module.
        """

        self.top_n = 5

        # find and load gather reports
        self.gather_raw_data = dict()
        for f in self.find_log_files("sourmash/gather", filehandles=True):
            d = gather2data(f)
            d.load_gather()
            if d.data:
                self.gather_raw_data[f["s_name"]] = d.data
            self.add_data_source(f, section="gather")

        self.gather_raw_data = self.ignore_samples(self.gather_raw_data)

        if len(self.gather_raw_data) == 0:
            return 0

        # data is in wrong format for writing to file
        # self.write_data_file(self.gather_raw_data, "gather")

        log.info("Found {} gather results".format(len(self.gather_raw_data)))

        # initialize variables to store summarized information
        self.gather_pct_per_match_all_samples = dict()
        self.gather_pct_unclassified_per_sample = dict()
        self.gather_pct_top_five_per_sample = dict()
        self.gather_top_five_matches = []

        # run functions to summarize information
        self.calculate_pct_per_match_all_samples()
        self.calculate_pct_unclassified_per_sample()
        # self.calculate_top_five_matches()
        self.calculate_pct_top_five_per_sample()

        # run functions to build multiqc report components
        self.general_stats_columns()
        self.top_five_barplot()

        # return the number of reports id'd, which will be logged in sourmash.py
        return len(self.gather_raw_data)

    def calculate_pct_per_match_all_samples(self):
        """
        sum the percent of each genome match across queries.
        e.g. if e. coli k12 is found to be 10% in sample A and 20% in sample B, the output would be e. coli k12 : 30.
        the output is a dictionary, gather_pc_per_match_all_samples, with match_names (genomes) as keys and summed percents as values.
        these values are used to identify the top 5 matches identified across all samples.
        """
        for s_name, data in self.gather_raw_data.items():
            for row in data:
                # convienence vars to make code easier to read
                match_name = row["match_name"]
                # loop over all samples and sum different variables
                if match_name not in self.gather_pct_per_match_all_samples:
                    self.gather_pct_per_match_all_samples[match_name] = 0
                self.gather_pct_per_match_all_samples[match_name] += row["pct_unique_weighted"]

    def calculate_pct_unclassified_per_sample(self):
        """
        calculate the percent of each sample that was unclassified.
        the output is a dictionary, gather_pct_unclassified_per_sample, with s_name (samples) as keys and the percent of the sample unclassified as values.
        """
        for s_name, data in self.gather_raw_data.items():
            for row in data:
                # sum over all matches for a given query to calculate the percent classified
                if s_name not in self.gather_pct_unclassified_per_sample:
                    self.gather_pct_unclassified_per_sample[s_name] = 0
                self.gather_pct_unclassified_per_sample[s_name] += row["pct_unique_weighted"]
            # convert to % unclassified
            self.gather_pct_unclassified_per_sample[s_name] = 100 - self.gather_pct_unclassified_per_sample[s_name]

    def calculate_pct_top_five_per_sample(self):
        """
        calculate the percent of each sample that is attributable to the top 5 genomes across all samples.
        the output is a dictionary, gather_pct_top_five_per_sample, with s_name (samples) as keys as the summed percent for the top 5 genomes as values.
        """
        # get top genomes matched across samples
        sorted_pct = sorted(self.gather_pct_per_match_all_samples.items(), key=lambda x: x[1], reverse=True)
        for pct_sum in sorted_pct[: self.top_n]:
            self.gather_top_five_matches.append(pct_sum[0])  # append the genome name to the top five list

        # calculate the pct attributable to the top 5 matches per sample
        for match_name in self.gather_top_five_matches:
            for s_name, data in self.gather_raw_data.items():
                if s_name not in self.gather_pct_top_five_per_sample:
                    self.gather_pct_top_five_per_sample[s_name] = 0
                for row in data:
                    if row["match_name"] == match_name:
                        self.gather_pct_top_five_per_sample[s_name] += row["pct_unique_weighted"]

    def general_stats_columns(self):
        """
        add columns to the general statistics table for % top 5 matches and % unclassified
        """
        # Column headers
        headers = OrderedDict()
        headers["% Top 5"] = {
            "title": "% Top 5 Genomes",
            "description": "Percentage of sample that was classified as one of the top 5 genomes ({})".format(
                ", ".join(self.gather_top_five_matches)
            ),
            "suffix": "%",
            "max": 100,
            "scale": "PuBu",
        }

        headers["% Unclassified"] = {
            "title": "% Unclassified",
            "description": "Percentage of sample that was unclassified",
            "suffix": "%",
            "max": 100,
            "scale": "OrRd",
        }

        # set table data
        tdata = {}
        for s_name, data in self.gather_raw_data.items():
            tdata[s_name] = {}
            tdata[s_name]["% Unclassified"] = self.gather_pct_unclassified_per_sample[s_name]
            tdata[s_name]["% Top 5"] = self.gather_pct_top_five_per_sample[s_name]

        self.general_stats_addcols(tdata, headers)

    def top_five_barplot(self):
        """Add a bar plot showing the percentage of top-5 genomes, the percentage of other genomes, and the unclassified percentage"""

        pd = []  # plot data
        cats = (
            list()
        )  # a list of categories to be shown as a color on the plot. includes the top 5 genomes, other, and unclassified
        pconfig = {
            "id": "gather-topfive-plot",
            "title": "sourmash: gather top genomes",
            "ylab": "percentages",
            "cpswitch": False,  # do not show the 'Counts / Percentages' switch, since gather only reports percentages
            "cpswitch_c_active": False,  # Initial display should show percentage, not counts
        }

        # create dictionaries for plot fill names and data (match percentages)
        match_cats = OrderedDict()
        match_data = dict()

        pct_shown = {}
        for match_name in self.gather_top_five_matches:
            # make match_name a fill category for the plot
            match_cats[match_name] = {"name": match_name}
            # pull out percents for this match_name from each sample
            for s_name, data in self.gather_raw_data.items():
                if s_name not in match_data:
                    match_data[s_name] = dict()
                if s_name not in pct_shown:
                    pct_shown[s_name] = 0
                for row in data:
                    if row["match_name"] == match_name:
                        if match_name not in match_data[s_name]:
                            match_data[s_name][match_name] = 0
                        match_data[s_name][match_name] += row["pct_unique_weighted"]
                        pct_shown[s_name] += row["pct_unique_weighted"]

        # Add in unclassified pct and other accounted for by genomes matches not in the top 5
        for s_name, data in self.gather_raw_data.items():
            # unclassified
            match_data[s_name]["unclassified"] = self.gather_pct_unclassified_per_sample[s_name]
            pct_shown[s_name] += self.gather_pct_unclassified_per_sample[s_name]
            # other
            match_data[s_name]["other"] = 100 - (
                self.gather_pct_unclassified_per_sample[s_name] + self.gather_pct_top_five_per_sample[s_name]
            )

        match_cats["other"] = {"name": "Other", "color": "#cccccc"}
        match_cats["unclassified"] = {"name": "Unclassified", "color": "#d4949c"}

        cats.append(match_cats)
        pd.append(match_data)

        self.add_section(
            name="Top gather genomes",
            anchor="gather-topfive",
            description=f"The percentage of the sample falling into the top {self.top_n} genome matches.",
            helptext=f"""
                To make this plot, the percentage of each sample assigned to a given genome is summed across all samples.
                The percentages for these top five genomes are then plotted, as well as the unclassified percentage.
                The unclassified count is always shown across all taxa ranks.
                The category _"Other"_ shows the difference between the above total percentages and the sum of the percentages
                in the top {self.top_n} genomes shown + unclassified. This should cover all genomes _not_ in the top {self.top_n}, +/- any rounding errors.
            """,
            plot=bargraph.plot(pd, cats, pconfig),
        )


class gather2data:
    """class to read in and parse the gather csv into a list of dictionaries."""

    def __init__(self, gather_file):
        self.data = []

        self.gather_file = gather_file

    def load_gather(self):
        gatherr = csv.DictReader(self.gather_file["f"])
        for line in gatherr:
            row = {
                "query_name": line["query_name"],
                "match_name": line["name"],
                "pct_unique_weighted": float(line["f_unique_weighted"]) * 100,
                "pct_match": float(line["f_match"]) * 100,
            }
            self.data.append(row)
