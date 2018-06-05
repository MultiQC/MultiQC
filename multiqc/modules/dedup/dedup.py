#!/usr/bin/env python

""" MultiQC module to parse output from DeDup """

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import re

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
        for f in self.find_log_files('dedup'):
            self.parse_dedup_log(f)

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


    def parse_dedup_log(self, f):
        regexes = {
            'input_file': r"Input file:\s+(\S+)",
            'total_reads': r"Total reads:\s+(\d+)",
            'reverse_removed': r"Reverse removed:\s+(\d+)",
            'forward_removed': r"Forward removed:\s+(\d+)",
            'merged_removed': r"Merged removed:\s+(\d+)",
            'total_removed': r"Total removed:\s+(\d+)",
            'duplication_rate': r"Duplication Rate:\s+([\d\.]+)",
        }
        parsed_data = dict()
        for k, r in regexes.items():
            r_search = re.search(r, f['f'], re.MULTILINE)
            if r_search:
                try:
                    parsed_data[k] = float(r_search.group(1))
                except ValueError:
                    parsed_data[k] = r_search.group(1)
        try:
            parsed_data['not_removed'] = parsed_data['total_reads'] - parsed_data['reverse_removed'] - parsed_data['forward_removed'] - parsed_data['merged_removed']
        except KeyError:
            log.debug('Could not calculate "not_removed"')

        if len(parsed_data) > 0:
            s_name = self.clean_s_name(os.path.basename(f['root']), f['root'])
            if 'input_file' in parsed_data:
                s_name = self.clean_s_name(parsed_data['input_file'], f['root'])
            self.dedup_data[s_name] = parsed_data

    def dedup_general_stats_table(self):
        """ Take the parsed stats from the DeDup report and add it to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['duplication_rate'] = {
            'title': 'Duplication Rate',
            'description': 'Percentage of reads categorised as a technical duplicate',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'OrRd',
            'format': '{:,.0f}',
            'modify': lambda x: x * 100.0
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
