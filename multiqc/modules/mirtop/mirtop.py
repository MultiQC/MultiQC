#!/usr/bin/env python

""" MultiQC module to parse output from mirtop"""

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import bargraph, beeswarm
from multiqc.modules.base_module import BaseMultiqcModule

import json

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='miRTop',
        anchor='mirtop', target='mirtop',
        href='https://github.com/miRTop',
        info="is a Command line tool to annotate miRNAs and isomiRs using a standard naming."
        )

        # Find and load any mirtop reports
        self.mirtop_data = dict()
        self.mirtop_keys = list()
        for c_file in self.find_log_files('mirtop'):
            content = json.loads(c_file['f'])
            self.parse_mirtop_report(content)
        # Filter to strip out ignored sample names
        self.mirtop_data = self.ignore_samples(self.mirtop_data)

        if len(self.mirtop_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.mirtop_data)))

        # Write parsed report data to a file
        self.write_data_file(self.mirtop_data, 'multiqc_mirtop')

        # Create very basic summary table
        self.mirtop_stats_table()

        # Create comprehensive beeswarm plots of all stats 
        self.add_section (
             name = 'IsoMiR Summary Statistics',
             anchor = 'mirtop-stats',
             description = "Detailed summary stats",
             plot = beeswarm.plot(self.mirtop_data)
         )



    def parse_mirtop_report (self, content):
        """ Parse the mirtop log file. """
        
        file_names = list()
        parsed_data = dict()
        for sample_name in content['metrics'].keys():
            log.info("Importing sample " + sample_name)
            parsed_data = content['metrics'][sample_name]
            parsed_data['read_count'] = parsed_data['isomiR_sum'] + parsed_data['ref_miRNA_sum']
            parsed_data['isomiR_perc'] = (parsed_data['isomiR_sum'] / parsed_data['read_count'])*100
            self.mirtop_data[sample_name] = parsed_data

    def mirtop_stats_table(self):
        """ Take the parsed stats from the mirtop report and add them to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['isomiR_sum'] = {
            'title': 'IsomiR reads',
            'description': 'read count summed over all isomiRs in sample',
            'scale': 'PuBu',
            'shared_key': 'read_co'
        }
        headers['ref_miRNA_sum'] = {
            'title': 'Reference reads',
            'description': 'read count summed over all reads mapping to the reference form of a miRNA',
            'scale': 'RdYlGn',
            'shared_key': 'read_co'
        }
        headers['read_count'] = {
            'title': 'Total reads',
            'description': 'all aligned reads',
            'scale': 'PuBu',
            'shared_key': 'read_co'
        } 
        headers['isomiR_perc'] = {
            'title': 'IsomiR %',
            'description': 'percentage of reads mapping to non-canonical forms of a microRNA',
            'min':0,
            'max':100,
            'suffix':'%',
            'scale': 'RdYlGn'
        }

        self.general_stats_addcols(self.mirtop_data, headers)
