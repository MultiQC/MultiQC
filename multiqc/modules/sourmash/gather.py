#!/usr/bin/env python

""" MultiQC module to parse similarity matrix output by sourmash gather """

import logging
import csv

from multiqc.plots import heatmap

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
        for f in self.find_log_files("sourmash/gather", filehandles=True)
            d = gather2data(f)
            if d.data:
                self.gather_raw_data[f["s_name"]] = d
            self.add_data_source(f, section="gather")

        self.gather_raw_data = self.ignore_samples(self.gather_raw_data)

        if len(self.gather_raw_data) == 0:
            raise UserWarning

        # data is in wrong format for writing to file
        # self.write_data_file(self.gather_raw_data, "gather")

        log.info("Found {} reports".format(len(self.gather_raw_data)))

        # initialize variables to store summarized information
        self.gather_pct_per_match_all_samples = dict()
        self.gather_pct_unclassified_per_sample = dict()
        self.gather_pct_top_five_per_sample = dict()
        self.gather_top_five_matches = []

        # run functions to summarize information
        self.calculate_pct_per_match_all_samples()
        self.calculate_pct_unclassified_per_sample()
        self.calculate_top_five_matches()
        self.calculate_pct_top_five_per_sample()

        # run functions to build multiqc report components
        self.general_stats_cols()
        #self.top_five_barplot()


    def calculate_pct_per_match_all_samples(self):
        """ 
        sum the percent of each genome match across queries.
        these values will be used to identify the top 5 matches.
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
        for s_name, data in self.gather_raw_data.items():
            for row in data:
            # sum over all matches for a given query to calculate the percent classified
                if s_name not in self.gather_pct_unclassified_per_sample:
                    self.gather_pct_unclassified_per_sample[s_name] = 0
                self.gather_pct_unclassified_per_sample[s_name] += row['pct_unique_weighted']
            # convert to % unclassified
            self.gather_pct_unclassified_per_sample[s_name] = 100 - self.gather_pct_unclassified_per_sample[s_name]

    def calculate_pct_top_five_per_sample(self):
    """
    calculate the percent of each sample that is attributable to the top 5 genomes across all samples
    """
        # get top genomes matched across samples
        sorted_pct = sorted(self.gather_pct_per_match_all_samples.items(), key=lambda x: x[1], reverse=True)
        for pct_sum in sorted_pct[: top_n]:
            self.gather_top_five_matches.append(pct_sum[0]) # append the genome name to the top five list

        # calculate the pct attributable to the top 5 matches per sample
        for match_name in self.gather_top_five_matches:
            for s_name, data in self.gather_raw_data.items():
                if s_name not in self.gather_pct_top_five_per_sample:
                    self.gather_pct_top_five_per_sample[s_name] = 0
                for row in data:
                    if row["match_name"] == match_name:
                        self.gather_pct_top_five_per_sample[s_name] += row['pct_unique_weighted']
   
    def general_stats_column(self):
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

class gather2data:
    """ class to read in and parse the gather csv into a list of dictionaries. """
    def __init__(self, gather_file):
        self.data = []

        self.load_gather(gather_file["f"])

    def load_gather(self, f):
        gatherr = csv.DictReader(f)
        for line in gatherr:    
            row = {
                "query_name" : line['query_name'],
                "match_name" : line['name'],
                "pct_unique_weighted" : float(line['f_unique_weighted']) * 100,
                "pct_match" : float(line['f_match']) * 100
            }
            self.data.append(row)
