#!/usr/bin/env python

""" MultiQC module to parse output from Bowtie 2 """

from __future__ import print_function
from collections import OrderedDict
import io
import json
import logging
import mmap
import os
import re

import multiqc

# Initialise the logger
log = logging.getLogger('MultiQC : {0:<14}'.format('Bowtie 2'))

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self, report):

        # Initialise the parent object
        super(MultiqcModule, self).__init__()

        # Static variables
        self.name = "Bowtie 2"
        self.anchor = "bowtie2"
        self.intro = '<p><a href="http://bowtie-bio.sourceforge.net/bowtie2/" target="_blank">Bowtie 2</a> \
            is ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences.</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']

        # Find and load any Bowtie 2 reports
        self.bowtie2_data = dict()
        for root, dirnames, filenames in os.walk(self.analysis_dir, followlinks=True):
            for fn in filenames:
                try:
                    if os.path.getsize(os.path.join(root,fn)) < 200000:
                        with io.open (os.path.join(root,fn), "r", encoding='utf-8') as f:
                            s = f.read()
                            parsed_data = self.parse_bowtie2_logs(s)
                            if parsed_data is not None:
                                s_name = self.clean_s_name(fn, root, prepend_dirs=report['prepend_dirs'])
                                if s_name in self.bowtie2_data:
                                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                                self.bowtie2_data[s_name] = parsed_data
                except (OSError, ValueError, UnicodeDecodeError):
                    log.debug("Couldn't read file when looking for output: {}".format(fn))

        if len(self.bowtie2_data) == 0:
            log.debug("Could not find any reports in {}".format(self.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.bowtie2_data)))

        # Write parsed report data to a file
        with io.open (os.path.join(self.output_dir, 'report_data', 'multiqc_bowtie2.txt'), "w", encoding='utf-8') as f:
            print( self.dict_to_csv( { k: { j: x for j, x in v.items() if j != 't_lengths'} for k, v in self.bowtie2_data.items() } ), file=f)

        self.sections = list()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.bowtie2_general_stats_table(report)

        # Alignment Rate Plot
        # Only one section, so add to the intro
        self.intro += self.bowtie2_alignment_plot()


    def parse_bowtie2_logs(self, s):
        # Check that this isn't actually Bismark using bowtie
        if s.find('Using bowtie 2 for aligning with bismark.', 0) >= 0: return None
        i = s.find('reads; of these:', 0)
        parsed_data = {}
        if i >= 0:
            regexes = {
                'reads_processed': r"(\d+) reads; of these:",
                'reads_aligned': r"(\d+) \([\d\.]+%\) aligned (?:concordantly )?exactly 1 time",
                'reads_aligned_percentage': r"\(([\d\.]+)%\) aligned (?:concordantly )?exactly 1 time",
                'not_aligned': r"(\d+) \([\d\.]+%\) aligned (?:concordantly )?0 times",
                'not_aligned_percentage': r"\(([\d\.]+)%\) aligned (?:concordantly )?0 times",
                'multimapped': r"(\d+) \([\d\.]+%\) aligned (?:concordantly )?>1 times",
                'multimapped_percentage': r"\(([\d\.]+)%\) aligned (?:concordantly )?>1 times",
                'overall_aligned_rate': r"([\d\.]+)% overall alignment rate",
            }

            for k, r in regexes.items():
                match = re.search(r, s)
                if match:
                    parsed_data[k] = float(match.group(1).replace(',', ''))
            
        if len(parsed_data) == 0: return None
        parsed_data['reads_other'] = parsed_data['reads_processed'] - parsed_data.get('reads_aligned', 0) - parsed_data.get('not_aligned', 0) - parsed_data.get('multimapped', 0)
        return parsed_data


    def bowtie2_general_stats_table(self, report):
        """ Take the parsed stats from the Bowtie 2 report and add it to the
        basic stats table at the top of the report """

        report['general_stats']['headers']['bowtie2_aligned'] = '<th class="chroma-col" data-chroma-scale="OrRd-rev" data-chroma-max="100" data-chroma-min="20"><span data-toggle="tooltip" title="Bowtie 2: overall alignment rate">% Aligned</span></th>'
        for samp, vals in self.bowtie2_data.items():
            report['general_stats']['rows'][samp]['bowtie2_aligned'] = '<td class="text-right">{:.1f}%</td>'.format(vals['overall_aligned_rate'])

    def bowtie2_alignment_plot (self):
        """ Make the HighCharts HTML to plot the alignment rates """
        
        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['reads_aligned'] = { 'color': '#8bbc21', 'name': '1 Alignment' }
        keys['multimapped'] =   { 'color': '#2f7ed8', 'name': '>1 Alignments' }
        keys['not_aligned'] =   { 'color': '#0d233a', 'name': 'Not aligned' }
        keys['reads_other'] =   { 'color': '#fd0000', 'name': 'Other' }
        
        # Config for the plot
        config = {
            'title': 'Bowtie 2 Alignment Scores',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        
        return self.plot_bargraph(self.bowtie2_data, keys, config)
