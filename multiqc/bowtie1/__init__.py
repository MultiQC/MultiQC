#!/usr/bin/env python

""" MultiQC module to parse output from Bowtie 1 """

from __future__ import print_function
from collections import OrderedDict
import json
import logging
import mmap
import os
import re

import multiqc
from multiqc import config

# Initialise the logger
log = logging.getLogger('MultiQC : {0:<14}'.format('Bowtie 1'))

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__()

        # Static variables
        self.name = "Bowtie"
        self.anchor = "bowtie1"
        self.intro = '<p><a href="http://bowtie-bio.sourceforge.net/" target="_blank">Bowtie 1</a> \
            is an ultrafast, memory-efficient short read aligner.</p>'

        # Find and load any Bowtie reports
        self.bowtie_data = dict()
        for f in self.find_log_files(contents_match='# reads processed:'):
            parsed_data = self.parse_bowtie_logs(f['f'])
            if parsed_data is not None:
                if f['s_name'] in self.bowtie_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.bowtie_data[f['s_name']] = parsed_data

        if len(self.bowtie_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.bowtie_data)))

        # Write parsed report data to a file
        self.write_csv_file(self.bowtie_data, 'multiqc_bowtie1.txt')

        self.sections = list()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.bowtie_general_stats_table()

        # Alignment Rate Plot
        # Only one section, so add to the intro
        self.intro += self.bowtie_alignment_plot()


    def parse_bowtie_logs(self, s):
        # Check that this isn't actually Bismark using bowtie
        if s.find('Using bowtie 1 for aligning with bismark.', 0) >= 0: return None
        parsed_data = {}
        regexes = {
            'reads_processed': r"# reads processed:\s+(\d+)",
            'reads_aligned': r"# reads with at least one reported alignment:\s+(\d+)",
            'reads_aligned_percentage': r"# reads with at least one reported alignment:\s+\d+\s+\(([\d\.]+)%\)",
            'not_aligned': r"# reads that failed to align:\s+(\d+)",
            'not_aligned_percentage': r"# reads that failed to align:\s+\d+\s+\(([\d\.]+)%\)",
            'multimapped': r"# reads with alignments suppressed due to -m:\s+(\d+)",
            'multimapped_percentage': r"# reads with alignments suppressed due to -m:\s+\d+\s+\(([\d\.]+)%\)"
        }
        for k, r in regexes.items():
            match = re.search(r, s)
            if match:
                parsed_data[k] = float(match.group(1).replace(',', ''))
        if len(parsed_data) == 0:
            return None
        return parsed_data


    def bowtie_general_stats_table(self):
        """ Take the parsed stats from the Bowtie report and add it to the
        basic stats table at the top of the report """

        config.general_stats['headers']['bowtie_aligned'] = '<th class="chroma-col" data-chroma-scale="OrRd-rev" data-chroma-max="100" data-chroma-min="20"><span data-toggle="tooltip" title="Bowtie 1: % reads with at least one reported alignment">% Aligned</span></th>'
        for samp, vals in self.bowtie_data.items():
            config.general_stats['rows'][samp]['bowtie_aligned'] = '<td class="text-right">{:.1f}%</td>'.format(vals['reads_aligned_percentage'])

    def bowtie_alignment_plot (self):
        """ Make the HighCharts HTML to plot the alignment rates """
        
        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['reads_aligned'] = { 'color': '#8bbc21', 'name': 'Aligned' }
        keys['multimapped'] =   { 'color': '#2f7ed8', 'name': 'Multimapped' }
        keys['not_aligned'] =   { 'color': '#0d233a', 'name': 'Not aligned' }
        
        # Config for the plot
        config = {
            'title': 'Bowtie 1 Alignment Scores',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        
        return self.plot_bargraph(self.bowtie_data, keys, config)
