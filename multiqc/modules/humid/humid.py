#!/usr/bin/env python

""" MultiQC module to parse output from Lima """

import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialse the parent object
        super(MultiqcModule, self).__init__(
            name="HUMID",
            anchor="humid",
            href="https://github.com/jfjlaros/dedup",
            info=" -- Error tolerant UMI aware FastQ deduplicator.",
            # No publication / DOI // doi=
        )

        # To store the summary data
        self.humid = dict()

        # Parse the output files
        self.parse_stat_files()

        # Remove filtered samples
        self.humid = self.ignore_samples(self.humid)

        # Let MultiQC know this module found no data
        if not self.humid: 
            raise UserWarning

        self.write_data_file(self.humid, "multiqc_humid")
        log.info(f"Found {len(self.humid)} reports")

    def parse_stat_files(self):
        for f in self.find_log_files("humid", filehandles=True):
            data = parse_stat_file(f["f"])
            # There is no sample name in the log, so we use the root of the
            # file as sample name (since the filename is always stats.dat
            name = f["root"]
            f["s_name"] = name
            self.humid[name] = data
            self.add_data_source(f)

def parse_stat_file(fin):
    """ Parse the stats file """
    data = dict()
    for line in fin:
        field, value = line.strip().split(': ')
        data[field] = int(value)
    process_stats(data)
    return data

def process_stats(stats):
    """ Process the statistics, to calculate some useful values """
    stats['filtered'] = stats['total'] - stats['usable']
    stats['duplicates'] = stats['total'] - stats['clusters'] - stats['filtered']
    # Sanity check
    assert stats['duplicates'] + stats['clusters'] + stats['filtered'] == stats['total']
