#!/usr/bin/env python

""" MultiQC module to parse output from Picard """

from __future__ import print_function
from collections import defaultdict, OrderedDict
import io
import json
import logging
import mmap
import os

import multiqc

# Initialise the logger
log = logging.getLogger('MultiQC : {0:<14}'.format('Picard'))

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self, report):

        # Initialise the parent object
        super(MultiqcModule, self).__init__()

        # Static variables
        self.name = "Picard"
        self.anchor = "picard"
        self.intro = '<p><a href="http://broadinstitute.github.io/picard/" target="_blank">Picard</a> \
            is a set of Java command line tools for manipulating high-throughput sequencing data.</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']

        # Find and load any Picard reports
        self.picard_dupMetrics_data = defaultdict(lambda: dict())
        for root, dirnames, filenames in os.walk(self.analysis_dir, followlinks=True):
            for fn in filenames:
                try:
                    if os.path.getsize(os.path.join(root,fn)) < 200000:
                        with io.open (os.path.join(root,fn), "r", encoding='utf-8') as f:
                            s = f.readlines()
                            for idx, l in enumerate(s):
                                if 'picard.sam.DuplicationMetrics' in l:
                                    s_name = fn
                                    s_name = s_name.split(".metrics",1)[0]
                                    s_name = self.clean_s_name(s_name, root, prepend_dirs=report['prepend_dirs'])
                                    if s_name in self.picard_dupMetrics_data:
                                        log.debug("Duplicate DuplicationMetrics sample name found! Overwriting: {}".format(s_name))
                                    keys = s[idx+1].split("\t")
                                    vals = s[idx+2].split("\t")
                                    for i, k in enumerate(keys):
                                        try:
                                            self.picard_dupMetrics_data[s_name][k] = float(vals[i])
                                        except ValueError:
                                            self.picard_dupMetrics_data[s_name][k] = vals[i]
                                    break # TODO: Deal with multiple libraries? Proper regex?
                except (OSError, ValueError, UnicodeDecodeError):
                    log.debug("Couldn't read file when looking for output: {}".format(fn))

        if len(self.picard_dupMetrics_data) == 0:
            log.debug("Could not find any reports in {}".format(self.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.picard_dupMetrics_data)))

        # Write parsed report data to a file
        with io.open (os.path.join(self.output_dir, 'report_data', 'multiqc_picard.txt'), "w", encoding='utf-8') as f:
            print( self.dict_to_csv( self.picard_dupMetrics_data ), file=f)

        self.sections = list()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.picard_stats_table(report)

        # Section 1 - Column chart of alignment stats
        self.sections.append({
            'name': 'Mark Duplicates',
            'anchor': 'picard-markduplicates',
            'content': self.mark_duplicates_plot()
        })


    def picard_stats_table(self, report):
        """ Take the parsed stats from the Picard report and add them to the
        basic stats table at the top of the report """

        report['general_stats']['headers']['picard_percent_duplication'] = '<th class="chroma-col" data-chroma-scale="OrRd" data-chroma-max="100" data-chroma-min="0"><span data-toggle="tooltip" title="Picard MarkDuplicates: Percent&nbsp;Duplication">%&nbsp;Dups</span></th>'
        for sn, data in self.picard_dupMetrics_data.items():
            report['general_stats']['rows'][sn]['picard_percent_duplication'] = '<td class="text-right">{:.1f}%</td>'.format(data['PERCENT_DUPLICATION']*100)

    def mark_duplicates_plot (self):
        """ Make the HighCharts HTML to plot the alignment rates """
        
        # NOTE: I had a hard time getting these numbers to add up as expected.
        # If you think I've done something wrong, let me know! Please add an
        # issue here: https://github.com/ewels/MultiQC/issues
        for sn in self.picard_dupMetrics_data.keys():
            self.picard_dupMetrics_data[sn]['UNPAIRED_READ_UNIQUE'] = self.picard_dupMetrics_data[sn]['UNPAIRED_READS_EXAMINED'] - self.picard_dupMetrics_data[sn]['UNPAIRED_READ_DUPLICATES']
            self.picard_dupMetrics_data[sn]['READ_PAIR_NOT_OPTICAL_DUPLICATES'] = self.picard_dupMetrics_data[sn]['READ_PAIR_DUPLICATES'] - self.picard_dupMetrics_data[sn]['READ_PAIR_OPTICAL_DUPLICATES']
            self.picard_dupMetrics_data[sn]['READ_PAIR_UNIQUE'] = self.picard_dupMetrics_data[sn]['READ_PAIRS_EXAMINED'] - self.picard_dupMetrics_data[sn]['READ_PAIR_DUPLICATES']
        
        keys = OrderedDict()
        keys_r = ['READ_PAIR_UNIQUE', 'UNPAIRED_READ_UNIQUE', 'READ_PAIR_NOT_OPTICAL_DUPLICATES',
                'READ_PAIR_OPTICAL_DUPLICATES', 'UNPAIRED_READ_DUPLICATES', 'UNMAPPED_READS']
        for k in keys_r:
            keys[k] = {'name': k.replace('_',' ').title()}
        
        # Config for the plot
        config = {
            'title': 'Picard Deduplication Stats',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        
        return self.plot_bargraph(self.picard_dupMetrics_data, keys, config)
