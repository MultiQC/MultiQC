#!/usr/bin/env python

""" MultiQC module to parse output from DeDup """

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import re
import json

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ DeDup module """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='DeDup', anchor='dedup',
        href="http://www.github.com/apeltzer/DeDup",
        info="is a tool for duplicate removal for merged/collapsed reads in ancient DNA analysis.")

        # Find and load any DeDup reports
        self.dedup_data = dict()

        for f in self.find_log_files('dedup',filehandles=True):
            self.parseJSON(f)

        # Filter to strip out ignored sample names
        self.dedup_data = self.ignore_samples(self.dedup_data)

        if len(self.dedup_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.dedup_data)))

        # Write parsed report data to a file
        self.write_data_file(self.dedup_data, 'multiqc_dedup')

        # Basic Stats Table
        self.dedup_general_stats_table()

        # Alignment Rate Plot
        self.add_section(
            description = 'This plot shows read categories that were either not removed (unique reads) or removed (duplicates).',
            plot = self.dedup_alignment_plot()
        )

    #Parse our nice little JSON file
    def parseJSON(self, f):
        """ Parse the JSON output from DeDup and save the summary statistics """
        try:
            parsed_json = json.load(f['f'])
        except JSONDecodeError as e:
            log.warn("Could not parse DeDup JSON: '{}'".format(f['fn']))
            print(e)
            return None

        #Get sample name from JSON first
        s_name = self.clean_s_name(parsed_json['metadata']['sample_name'],'')
        self.add_data_source(f, s_name)

        metrics_dict = parsed_json['metrics']

        for k in metrics_dict:
            metrics_dict[k] = float(metrics_dict[k])

        #Compute not removed from given values
        metrics_dict['not_removed'] = metrics_dict['total_reads'] - metrics_dict['reverse_removed'] - metrics_dict['forward_removed'] - metrics_dict['merged_removed']

        #Add all in the main data_table
        self.dedup_data[s_name] = metrics_dict


    def dedup_general_stats_table(self):
        """ Take the parsed stats from the DeDup report and add it to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['dup_rate'] = {
            'title': 'Duplication Rate',
            'description': 'Percentage of reads categorised as a technical duplicate',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'OrRd',
            'format': '{:,.0f}',
            'modify': lambda x: x * 100.0
        }
        headers['clusterfactor'] = {
            'title': 'ClusterFactor',
            'description': 'CF~1 means high library complexity. Large CF means not worth sequencing deeper.',
            'min': 1,
            'max': 100,
            'scale': 'OrRd',
            'format': '{:,.2f}',
        }
        self.general_stats_addcols(self.dedup_data, headers)

    def dedup_alignment_plot (self):
        """ Make the HighCharts HTML to plot the duplication rates """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['not_removed'] = { 'name': 'Not Removed' }
        keys['reverse_removed'] = { 'name': 'Reverse Removed' }
        keys['forward_removed'] =   { 'name': 'Forward Removed' }
        keys['merged_removed'] =   { 'name': 'Merged Removed' }

        # Config for the plot
        config = {
            'id': 'dedup_rates',
            'title': 'DeDup: Deduplicated Reads',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads',
            'hide_zero_cats': False
        }

        return bargraph.plot(self.dedup_data, keys, config)
