#!/usr/bin/env python

""" MultiQC module to parse output from Picard """

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import re

from multiqc import config, BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Picard', anchor='picard', 
        href='http://broadinstitute.github.io/picard/', 
        info="is a set of Java command line tools for manipulating high-"\
        "throughput sequencing data.")

        # Find and load any Picard reports
        self.picard_dupMetrics_data = dict()
        for f in self.find_log_files(contents_match='picard.sam.MarkDuplicates', filehandles=True):
            self.parse_picard_dupMetrics(f)   

        if len(self.picard_dupMetrics_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.picard_dupMetrics_data)))

        # Write parsed report data to a file
        self.write_csv_file(self.picard_dupMetrics_data, 'multiqc_picard.txt')

        self.sections = list()

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.picard_stats_table()

        # Section 1 - Column chart of alignment stats
        self.sections.append({
            'name': 'Mark Duplicates',
            'anchor': 'picard-markduplicates',
            'content': self.mark_duplicates_plot()
        })


    def parse_picard_dupMetrics(self, f):
        """ Go through log file looking for picard output """
        s_name = None
        for l in f['f']:
            # New log starting
            if 'picard.sam.MarkDuplicates' in l and 'INPUT' in l:
                s_name = None
                
                # Pull sample name from input
                fn_search = re.search("INPUT=\[([^\]]+)\]", l)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1))
                    s_name = self.clean_s_name(s_name, f['root'])
            
            if s_name is not None:
                if 'picard.sam.DuplicationMetrics' in l and '## METRICS CLASS' in l:
                    if s_name in self.picard_dupMetrics_data:
                        log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
                    self.picard_dupMetrics_data[s_name] = dict()
                    keys = f['f'].readline().split("\t")
                    vals = f['f'].readline().split("\t")
                    for i, k in enumerate(keys):
                        try:
                            self.picard_dupMetrics_data[s_name][k] = float(vals[i])
                        except ValueError:
                            self.picard_dupMetrics_data[s_name][k] = vals[i]
                    s_name = None
        
        for s_name in self.picard_dupMetrics_data.keys():
            if len(self.picard_dupMetrics_data[s_name]) == 0:
                self.picard_dupMetrics_data.pop(s_name, None)
                log.debug("Removing {} as no data parsed".format(s_name))
        
    
    def picard_stats_table(self):
        """ Take the parsed stats from the Picard report and add them to the
        basic stats table at the top of the report """
        
        headers = OrderedDict()
        headers['PERCENT_DUPLICATION'] = {
            'title': '% Dups',
            'description': 'MarkDuplicates - Percent Duplication',
            'max': 100,
            'min': 0,
            'scale': 'OrRd',
            'format': '{:.1f}%',
            'modify': lambda x: float(x) * 100
        }
        self.general_stats_addcols(self.picard_dupMetrics_data, headers)


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
            'cpswitch_counts_label': 'Number of Reads',
            'cpswitch_c_active': False
        }
        
        return self.plot_bargraph(self.picard_dupMetrics_data, keys, config)
