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

        # sum counts across all samples, so that we can pick top 5
        self.gather_total_pct_unique_weighted = dict()
        #self.sum_sample_counts()

        #self.general_stats_cols()
        #self.top_five_barplot()


    def calculate_total_pct_unique(self):
        """ 
        sum the match fractions across queries.
        these values will be used to identify the top 5 matches with or without abundance weighting 
        """
        for s_name, data in self.gather_raw_data.items():
            for row in data:
                # convienence vars to make code easier to read
                match_name = row["match_name"] 
                # loop over all samples and sum different variables 
                if match_name not in self.gather_total_pct_unique_weighted:
                    self.gather_total_pct_unique_weighted[match_name] = 0
                self.gather_total_pct_unique_weighted[match_name] += row["pct_unique_weighted"]


    def general_stats_column(self):
        """
        add columns to the general statistics table for % top 5 matches and % unclassified
        """

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
