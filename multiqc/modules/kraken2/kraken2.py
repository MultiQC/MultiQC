#!/usr/bin/env python

""" MultiQC module to parse output from KRAKEN2 """

from __future__ import print_function
from collections import OrderedDict
import os
import logging
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ Kraken2 module """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='kraken2', anchor='kraken2',
        href="https://ccb.jhu.edu/software/kraken2/",
        info="Kraken 2 is the newest version of Kraken, a taxonomic classification system using exact k-mer matches to achieve high accuracy and fast classification speeds. This classifier matches each k-mer within a query sequence to the lowest common ancestor (LCA) of all genomes containing the given k-mer. The k-mer assignments inform the classification algorithm. ")

        # Find and load any kraken2 reports
        self.kraken2_data = dict()
        for f in self.find_log_files('kraken2'):
            self.parse_kraken2_log(f)

        # Filter to strip out ignored sample names
        self.kraken2_data = self.ignore_samples(self.kraken2_data)

        if len(self.kraken2_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.kraken2_data)))

        # Write parsed report data to a file
        self.write_data_file(self.kraken2_data, 'multiqc_kraken2')

        # Basic Stats Table
        self.kraken2_general_stats_table()


    def parse_kraken_log(self, f):
        kraken_data = OrderedDict()
        regex = "^\s{1,2}(\d{1,2}\.\d{1,2})\t(\d+)\t(\d+)\t([UDKPCOFGS-])\t(\d+)\s+(.+)"
        for l in f['f']:
            #TODO the regex parses 6 groups per matched line (
            # Example line:  11.66	98148	98148	U	0	unclassified
            # G1: 11.66, G2: 98148, G3: 98148, G4: U, G5: 0 , G6:Unclassified
            # Works on all 5 input files provided at MultiQC Testdata - will now parse this to a Dict and then think about vis
            # Worst case: Top Level info in main table and/or table info
        # Search regexes for stats
        for l in f['f'].splitlines():
            match = re.search(regex, l)
            if match:
                percent = float(match.group(1))
                counts = int(match.group(2))
                counts_2 = int(match.group(3))
                tax = match.group(4)
                counts_3 = int(match.group(5))
                classif = match.group(6)
                
                primers[primer] = counts
        log.info("Found {} primers".format(len(primers)))
        return primers


            

    