#!/usr/bin/env python

""" MultiQC module to parse output from Picard """

from __future__ import print_function
from collections import defaultdict, OrderedDict
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
        
        sp = config.sp['picard']

        #### MarkDuplicates reports
        self.picard_dupMetrics_data = dict()
        for f in self.find_log_files(sp['markdups'], filehandles=True):
            self.parse_picard_dupMetrics(f)
        
        #### CollectInsertSizeMetrics reports
        self.picard_insertSize_data = dict()
        self.picard_insertSize_histogram = dict()
        self.picard_insertSize_medians = dict()
        for f in self.find_log_files(sp['insertsize'], filehandles=True):
            self.parse_picard_insertSize(f)
        # Calculate summed median values for all read orientations
        for s_name in self.picard_insertSize_histogram:
            median = None
            j = 0
            for idx, c in self.picard_insertSize_histogram[s_name].items():
                j += c
                if j > (self.picard_insertSize_medians[s_name]['total_count'] / 2):
                    self.picard_insertSize_medians[s_name]['summed_median'] = idx
                    break
        
        #### CollectGcBias reports
        self.picard_GCbias_data = dict()
        for f in self.find_log_files(sp['gcbias'], filehandles=True):
            self.parse_picard_GCbiasMetrics(f)
        
        #### CalculateHsMetric
        self.picard_HsMetrics_data = dict()
        for f in self.find_log_files(sp['hsmetrics'], filehandles=True):
            self.parse_picard_HsMetrics(f)
        
        num_reports = ( len(self.picard_dupMetrics_data) + len(self.picard_insertSize_data) +
                len(self.picard_GCbias_data) + len(self.picard_HsMetrics_data) )
        
        if num_reports == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        self.sections = list()
        
        self.picard_stats_table()
        
        # Mark Duplicates data
        if len(self.picard_dupMetrics_data) > 0:
            log.info("Found {} dupMetrics reports".format(len(self.picard_dupMetrics_data)))
            self.write_csv_file(self.picard_dupMetrics_data, 'multiqc_picard_dups.txt')
            self.sections.append({
                'name': 'Mark Duplicates',
                'anchor': 'picard-markduplicates',
                'content': self.mark_duplicates_plot()
            })
        
        # Insert Size data
        if len(self.picard_insertSize_data) > 0:
            log.info("Found {} insertSize reports".format(len(self.picard_insertSize_data)))
            self.write_csv_file(self.picard_insertSize_data, 'multiqc_picard_insertSize.txt')
            self.sections.append({
                'name': 'Insert Size',
                'anchor': 'picard-insertsize',
                'content': '<p>Plot shows the number of reads at a given insert size. Reads with different orientations are summed.</p>' + 
                                self.insert_size_plot()
            })
        
        # GC Bias data
        if len(self.picard_GCbias_data) > 0:
            log.info("Found {} GCbias reports".format(len(self.picard_GCbias_data)))
            self.sections.append({
                'name': 'GC Coverage Bias',
                'anchor': 'picard-gcbias',
                'content': '<p>This plot shows bias in coverage across regions of the genome with varying GC content.'\
                    ' A perfect library would be a flat line at <code>y = 1</code>.</p>'+self.GCbias_plot()
            })
        
        # HsMetrics data
        if len(self.picard_HsMetrics_data) > 0:
            log.info("Found {} HsMetrics reports".format(len(self.picard_HsMetrics_data)))
            self.write_csv_file(self.picard_HsMetrics_data, 'multiqc_picard_HsMetrics.txt')


    def parse_picard_dupMetrics(self, f):
        """ Parse MarkDuplicate Picard Output """
        s_name = None
        for l in f['f']:
            # New log starting
            if 'picard.sam.MarkDuplicates' in l and 'INPUT' in l:
                s_name = None
                
                # Pull sample name from input
                fn_search = re.search("INPUT=\[?([^\\s]+)\]?", l)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1))
                    s_name = self.clean_s_name(s_name, f['root'])
            
            if s_name is not None:
                if 'picard.sam.DuplicationMetrics' in l and '## METRICS CLASS' in l:
                    if s_name in self.picard_dupMetrics_data:
                        log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
                    self.add_data_source(f, s_name, section='DuplicationMetrics')
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
                    ins = int(sections[0])
                    tot_count = sum( [int(x) for x in sections[1:]] )
                    self.picard_insertSize_histogram[s_name][ins] = tot_count
                    self.picard_insertSize_medians[s_name]['total_count'] += tot_count
                except ValueError:
                    # Reset in case we have more in this log file
                    s_name = None
                    in_hist = False
                    
            # New log starting
            if 'InsertSizeMetrics' in l and 'INPUT' in l:
                s_name = None
                # Pull sample name from input
                fn_search = re.search("INPUT=\[?([^\\s]+)\]?", l)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1))
                    s_name = self.clean_s_name(s_name, f['root'])
            
            if s_name is not None:
                if 'InsertSizeMetrics' in l and '## METRICS CLASS' in l:
                    if s_name in self.picard_insertSize_data:
                        log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
                    self.add_data_source(f, s_name, section='InsertSizeMetrics')
                    keys = f['f'].readline().split("\t")
                    vals = f['f'].readline().split("\t")
                    self.picard_insertSize_medians[s_name] = {'total_count': 0}
                    while len(vals) == len(keys):
                        pair_orientation = vals[7]
                        rowkey = '{}_{}'.format(s_name, pair_orientation)
                        self.picard_insertSize_data[rowkey] = OrderedDict()
                        self.picard_insertSize_data[rowkey]['SAMPLE_NAME'] = s_name
                        for i, k in enumerate(keys):
                            try:
                                self.picard_insertSize_data[rowkey][k] = float(vals[i])
                            except ValueError:
                                self.picard_insertSize_data[rowkey][k] = vals[i]
                        vals = f['f'].readline().split("\t")
                    
                    # Skip lines on to histogram
                    l = f['f'].readline()
                    l = f['f'].readline()
                    
                    self.picard_insertSize_histogram[s_name] = dict()
                    in_hist = True        
        
        for key in self.picard_insertSize_data.keys():
            if len(self.picard_insertSize_data[key]) == 0:
                self.picard_insertSize_data.pop(key, None)
        for s_name in self.picard_insertSize_histogram.keys():
            if len(self.picard_insertSize_histogram[s_name]) == 0:
                self.picard_insertSize_histogram.pop(s_name, None)
                self.picard_insertSize_medians.pop(s_name, None)
                log.debug("Removing {} as no data parsed".format(s_name))
        
    
    def parse_picard_GCbiasMetrics(self, f):
        """ Parse GCBiasMetrics Picard Output """
        s_name = None
        gc_col = None
        cov_col = None
        for l in f['f']:
            # New log starting
            if 'picard.analysis.CollectGcBiasMetrics' in l and 'INPUT' in l:
                s_name = None
                
                # Pull sample name from input
                fn_search = re.search("INPUT=\[?([^\\s]+)\]?", l)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1))
                    s_name = self.clean_s_name(s_name, f['root'])
            
            if s_name is not None:
                if gc_col is not None and cov_col is not None :
                    try:
                        # Note that GC isn't always the first column.
                        s = l.split("\t")
                        self.picard_GCbias_data[s_name][ int(s[gc_col]) ] = float(s[cov_col])
                    except IndexError:
                        s_name = None
                        gc_col = None
                        cov_col = None
                
                if 'picard.analysis.GcBiasDetailMetrics' in l and '## METRICS CLASS' in l:
                    if s_name in self.picard_GCbias_data:
                        log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
                    self.add_data_source(f, s_name, section='GcBiasDetailMetrics')
                    self.picard_GCbias_data[s_name] = dict()
                    # Get header - find columns with the data we want
                    l = f['f'].readline()
                    s = l.split("\t")
                    gc_col = s.index('GC')
                    cov_col = s.index('NORMALIZED_COVERAGE')
                    
        
        for s_name in self.picard_GCbias_data.keys():
            if len(self.picard_GCbias_data[s_name]) == 0:
                self.picard_GCbias_data.pop(s_name, None)
                log.debug("Removing {} as no data parsed".format(s_name))
    
    
    def parse_picard_HsMetrics(self, f):
        """ Parse HsMetric Picard Output """
        parsed_data = dict()
        s_name = None
        keys = None
        for l in f['f']:
            # New log starting
            if 'picard.analysis.directed.CalculateHsMetrics' in l and 'INPUT' in l:
                s_name = None
                keys = None
                
                # Pull sample name from input
                fn_search = re.search("INPUT=\[?([^\\s]+)\]?", l)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1))
                    s_name = self.clean_s_name(s_name, f['root'])
                    parsed_data[s_name] = dict()
            
            if s_name is not None:
                if 'picard.analysis.directed.HsMetrics' in l and '## METRICS CLASS' in l:
                    keys = f['f'].readline().split("\t")
                elif keys:
                    vals = l.split("\t")
                    if len(vals) == len(keys):
                        j = 'NA'
                        if keys[0] == 'BAIT_SET':
                            j = vals[0]
                        parsed_data[s_name][j] = dict()
                        for i, k in enumerate(keys):
                            try:
                                parsed_data[s_name][j][k] = float(vals[i])
                            except ValueError:
                                parsed_data[s_name][j][k] = vals[i]
                    else:
                        s_name = None
                        keys = None
        
        # Remove empty dictionaries
        for s_name in parsed_data.keys():
            for j in parsed_data[s_name].keys():
                if len(parsed_data[s_name][j]) == 0:
                    parsed_data[s_name].pop(j, None)
            if len(parsed_data[s_name]) == 0:
                parsed_data.pop(s_name, None)
        
        # Manipulate sample names if multiple baits found
        for s_name in parsed_data.keys():
            for j in parsed_data[s_name].keys():
                this_s_name = s_name
                if(len(parsed_data[s_name]) > 1):
                    this_s_name = "{}: {}".format(s_name, j)
                if this_s_name in self.picard_HsMetrics_data:
                    log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], this_s_name))
                self.add_data_source(f, this_s_name, section='HsMetrics')
                self.picard_HsMetrics_data[this_s_name] = parsed_data[s_name][j]
    
    
    def picard_stats_table(self):
        """ Take the parsed stats from Picard and add them to the
        basic stats table at the top of the report """
        
        data = defaultdict(lambda:dict())
        for s_name in self.picard_dupMetrics_data:
            data[s_name]['PERCENT_DUPLICATION'] = self.picard_dupMetrics_data[s_name]['PERCENT_DUPLICATION']
        for s_name in self.picard_insertSize_medians:
            data[s_name]['summed_median'] = self.picard_insertSize_medians[s_name]['summed_median']
        for s_name in self.picard_HsMetrics_data:
            data[s_name]['FOLD_ENRICHMENT'] = self.picard_HsMetrics_data[s_name]['FOLD_ENRICHMENT']
            if data[s_name]['FOLD_ENRICHMENT'] == '?':
                data[s_name]['FOLD_ENRICHMENT'] = -1
            data[s_name]['PCT_TARGET_BASES_30X'] = self.picard_HsMetrics_data[s_name]['PCT_TARGET_BASES_30X']
        
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
        headers['summed_median'] = {
            'title': 'Insert Size',
            'description': 'Median Insert Size, all read orientations (bp)',
            'min': 0,
            'format': '{:.0f}',
            'scale': 'GnBu',
        }
        headers['FOLD_ENRICHMENT'] = {
            'title': 'Fold Enrichment',
            'min': 0,
            'format': '{:.0f}',
            'scale': 'Blues',
        }
        headers['PCT_TARGET_BASES_30X'] = {
            'title': 'Target Bases 30X',
            'description': 'Percent of target bases with coverage &ge; 30X',
            'max': 100,
            'min': 0,
            'format': '{:.0f}%',
            'scale': 'RdYlGn',
            'modify': lambda x: float(x) * 100
        }
        self.general_stats_addcols(data, headers)


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

    def GCbias_plot(self):
        """ Plot the Picard GC Bias plot """
        
        pconfig = {
            'id': 'picard_gcbias_plot',
            'title': 'GC Coverage Bias',
            'ylab': 'Normalized Coverage',
            'xlab': '% GC',
            'xmin': 0,
            'xmax': 100,
            'xDecimals': False,
            'ymin': 0,
            'yCeiling': 10,
            'tt_label': '<b>{point.x} %GC</b>: {point.y:.2f}',
            'yPlotLines': [
                {'value': 1, 'color': '#999999', 'width': 2, 'dashStyle': 'LongDash'},
            ]
        }
        
        return self.plot_xy_data(self.picard_GCbias_data, pconfig)
    