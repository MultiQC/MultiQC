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
    self.picard_GCbias_data = dict()
    
    # Go through logs and find Metrics
    for f in self.find_log_files(config.sp['picard']['gcbias'], filehandles=True):
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
                        s = l.strip("\n").split("\t")
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
                    s = l.strip("\n").split("\t")
                    gc_col = s.index('GC')
                    cov_col = s.index('NORMALIZED_COVERAGE')
                    
        
        for s_name in self.picard_GCbias_data.keys():
            if len(self.picard_GCbias_data[s_name]) == 0:
                self.picard_GCbias_data.pop(s_name, None)
                log.debug("Removing {} as no data parsed".format(s_name))
    

    
    
    if len(self.picard_GCbias_data) > 0:        
        
        # Plot the graph
        
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
        self.sections.append({
            'name': 'GC Coverage Bias',
            'anchor': 'picard-gcbias',
            'content': '<p>This plot shows bias in coverage across regions of the genome with varying GC content.'\
                ' A perfect library would be a flat line at <code>y = 1</code>.</p>' + 
                plots.linegraph.plot(self.picard_GCbias_data, pconfig)
        })
    
    
    # Return the number of detected samples to the parent module
    return len(self.picard_GCbias_data)
    