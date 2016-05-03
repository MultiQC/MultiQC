#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard InsertSizeMetrics """

from collections import OrderedDict
import logging
import os
import re

from multiqc import config, plots

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find Picard InsertSizeMetrics reports and parse their data """
    
    # Set up vars
    self.picard_insertSize_data = dict()
    self.picard_insertSize_histogram = dict()
    self.picard_insertSize_medians = dict()
    
    # Go through logs and find Metrics
    for f in self.find_log_files(config.sp['picard']['insertsize'], filehandles=True):
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
                    keys = f['f'].readline().strip("\n").split("\t")
                    vals = f['f'].readline().strip("\n").split("\t")
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
                        vals = f['f'].readline().strip("\n").split("\t")
                    
                    # Skip lines on to histogram
                    l = f['f'].readline().strip("\n")
                    l = f['f'].readline().strip("\n")
                    
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
    
    # Calculate summed median values for all read orientations
    for s_name in self.picard_insertSize_histogram:
        median = None
        j = 0
        for idx, c in self.picard_insertSize_histogram[s_name].items():
            j += c
            if j > (self.picard_insertSize_medians[s_name]['total_count'] / 2):
                self.picard_insertSize_medians[s_name]['summed_median'] = idx
                break
    
    
    if len(self.picard_insertSize_data) > 0:
        
        # Write parsed data to a file
        self.write_data_file(self.picard_insertSize_data, 'multiqc_picard_insertSize')
        
        # Add to general stats table
        self.general_stats_headers['summed_median'] = {
            'title': 'Insert Size',
            'description': 'Median Insert Size, all read orientations (bp)',
            'min': 0,
            'suffix': 'bp',
            'format': '{:.0f}',
            'scale': 'GnBu',
        }
        for s_name in self.picard_insertSize_medians:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update( self.picard_insertSize_medians[s_name] )
        
        
        # Make a normalised percentage version of the data
        data_percent = {}
        for s_name, data in self.picard_insertSize_histogram.items():
            data_percent[s_name] = OrderedDict()
            total = float( sum( data.values() ) )
            for k, v in data.items():
                data_percent[s_name][k] = (v/total)*100
        
        # Plot the data and add section
        pconfig = {
            'smooth_points': 500,
            'smooth_points_sumcounts': [True, False],
            'id': 'picard_insert_size',
            'title': 'Insert Size',
            'ylab': 'Count',
            'xlab': 'Insert Size (bp)',
            'xDecimals': False,
            'tt_label': '<b>{point.x} bp</b>: {point.y:.0f}',
            'ymin': 0,
            'data_labels': [
                {'name': 'Counts', 'ylab': 'Coverage'},
                {'name': 'Percentages', 'ylab': 'Percentage of Counts'}
            ]
        }
        self.sections.append({
            'name': 'Insert Size',
            'anchor': 'picard-insertsize',
            'content': '<p>Plot shows the number of reads at a given insert size. Reads with different orientations are summed.</p>' + 
                        plots.linegraph.plot([self.picard_insertSize_histogram, data_percent], pconfig)
        })
    
    
    # Return the number of detected samples to the parent module
    return len(self.picard_insertSize_data)
    