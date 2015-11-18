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

        # Find and load any Picard MarkDuplicates reports
        self.picard_dupMetrics_data = dict()
        for f in self.find_log_files(contents_match='picard.sam.MarkDuplicates', filehandles=True):
            self.parse_picard_dupMetrics(f)
        
        # Find and load any Picard CollectInsertSizeMetrics reports
        self.picard_insertSize_data = dict()
        self.picard_insertSize_histogram = dict()
        for f in self.find_log_files(contents_match='picard.analysis.CollectInsertSizeMetrics', filehandles=True):
            self.parse_picard_insertSize(f)
        
        num_reports = len(self.picard_dupMetrics_data) + len(self.picard_insertSize_data)
        
        if num_reports == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(num_reports))
        
        self.sections = list()
        
        # Mark Duplicates data
        if len(self.picard_dupMetrics_data) > 0:
            self.write_csv_file(self.picard_dupMetrics_data, 'multiqc_picard_dups.txt')
            self.picard_stats_table_markdups()
            self.sections.append({
                'name': 'Mark Duplicates',
                'anchor': 'picard-markduplicates',
                'content': self.mark_duplicates_plot()
            })
        
        # Insert Size data
        if len(self.picard_insertSize_data) > 0:
            self.write_csv_file(self.picard_insertSize_data, 'multiqc_picard_insertSize.txt')
            self.picard_stats_table_insertSize()
            self.sections.append({
                'name': 'Insert Size',
                'anchor': 'picard-insertsize',
                'content': self.insert_size_plot()
            })


    def parse_picard_dupMetrics(self, f):
        """ Parse MarkDuplicate Picard Output """
        s_name = None
        for l in f['f']:
            # New log starting
            if 'picard.sam.MarkDuplicates' in l and 'INPUT' in l:
                s_name = None
                
                # Pull sample name from input
                fn_search = re.search("INPUT=\[?([^\\\s]+)\]?", l)
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
    
    
    
    def parse_picard_insertSize(self, f):
        """ Parse CollectInsertSizeMetrics Picard Output """
        s_name = None
        in_hist = False
        for l in f['f']:
            
            # Catch the histogram values
            if s_name is not None and in_hist is True:
                try:
                    sections = l.split("\t")
                    self.picard_insertSize_histogram[s_name][int(sections[0])] = int(sections[1])
                except IndexError:
                    s_name = None
                    in_hist = False
            
            # New log starting
            if 'InsertSizeMetrics' in l and 'INPUT' in l:
                s_name = None
                # Pull sample name from input
                fn_search = re.search("INPUT=\[?([^\\\s]+)\]?", l)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1))
                    s_name = self.clean_s_name(s_name, f['root'])
            
            if s_name is not None:
                if 'InsertSizeMetrics' in l and '## METRICS CLASS' in l:
                    if s_name in self.picard_insertSize_data:
                        log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
                    self.picard_insertSize_data[s_name] = dict()
                    keys = f['f'].readline().split("\t")
                    vals = f['f'].readline().split("\t")
                    for i, k in enumerate(keys):
                        try:
                            self.picard_insertSize_data[s_name][k] = float(vals[i])
                        except ValueError:
                            self.picard_insertSize_data[s_name][k] = vals[i]
                    
                    # Skip lines on to histogram
                    l = f['f'].readline()
                    l = f['f'].readline()
                    l = f['f'].readline()
                    
                    self.picard_insertSize_histogram[s_name] = dict()
                    in_hist = True
        
        for s_name in self.picard_insertSize_data.keys():
            if len(self.picard_insertSize_data[s_name]) == 0:
                self.picard_insertSize_data.pop(s_name, None)
                self.picard_insertSize_histogram.pop(s_name, None)
                log.debug("Removing {} as no data parsed".format(s_name))
        
    
    def picard_stats_table_markdups(self):
        """ Take the parsed stats from the Picard Mark Duplicates report and add them to the
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
    
    
    def picard_stats_table_insertSize(self):
        """ Take the parsed stats from the Picard Insert Size report and add them to the
        basic stats table at the top of the report """
        
        headers = OrderedDict()
        headers['MEDIAN_INSERT_SIZE'] = {
            'title': 'Insert Size',
            'description': 'Median Insert Size',
            'min': 0,
            'scale': 'GnBu',
        }
        self.general_stats_addcols(self.picard_insertSize_data, headers)


    def mark_duplicates_plot (self):
        """ Make the Picard Mark Duplicates scores """
        
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
    
    
    def insert_size_plot(self):
        """ Plot the Picard Insert Size histograms """
        
        pconfig = {
            'id': 'picard_insert_size',
            'title': 'Insert Size',
            'ylab': 'Count',
            'xlab': 'Insert Size (bp)',
            'xDecimals': False,
            'tt_label': '<b>{point.x} bp</b>: {point.y:.0f}',
            'ymin': 0,
        }
        
        return self.plot_xy_data(self.picard_insertSize_histogram, pconfig)
