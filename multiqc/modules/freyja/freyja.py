#!/usr/bin/env python

""" MultiQC module to parse output from Lima """

import logging
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Freyja',
            anchor='freyja',
            href="https://github.com/andersen-lab/Freyja",
            info="Recover relative lineage abundances from mixed SARS-CoV-2 samples."
        )

        # To store the summary data
        self.freyja_data = dict()

        # Parse the output files
        self.parse_stat_files()

        # Remove filtered samples
        self.freyja_data = self.ignore_samples(self.freyja_data)

        # Let MultiQC know this module found no data
        if len(self.freyja_data) == 0:
            raise UserWarning

        log.info(f"Found {len(self.freyja_data)} reports")
        self.write_data_file(self.freyja_data, "multiqc_freyja")

    def parse_stat_files(self, f):    
        for f in self.find_log_files('freyja'):
            s_name = self.clean_s_name(f["root"], f)
            data = parse_stat_file(f["f"], s_name)
            if data:
                # There is no sample name in the log, so we use the root of the
                # file as sample name (since the filename is always stats.dat
                if s_name in self.freyja_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.freyja_data[s_name] = data
                self.add_data_source(f, s_name)
                
# TODO: This needs to be customised for Freyja        
def parse_stat_file(fin, s_name):
    """Parse the stats file"""
    data = dict()
    for line in fin:
        field, value = line.strip().split(": ")
        data[field] = int(value)
    if process_stats(data, s_name):
        return data

# TODO: This needs to be customised for Freyja    
def process_stats(stats, s_name):
    """Process the statistics, to calculate some useful values"""
    stats["filtered"] = stats["total"] - stats["usable"]
    stats["duplicates"] = stats["total"] - stats["clusters"] - stats["filtered"]
    # Sanity check
    try:
        assert stats["duplicates"] + stats["clusters"] + stats["filtered"] == stats["total"]
    except AssertionError:
        log.warning(f"HUMID stats looked wrong, skipping: {s_name}")
        return False