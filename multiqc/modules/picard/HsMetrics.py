#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard HsMetrics """

from collections import OrderedDict
import logging
import os
import re

from multiqc import config
from multiqc import plots

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find Picard HsMetrics reports and parse their data """
    
    # Set up vars
    self.picard_HsMetrics_data = dict()
    
    # Go through logs and find Metrics
    for f in self.find_log_files(config.sp['picard']['hsmetrics'], filehandles=True):
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
                    keys = f['f'].readline().strip("\n").split("\t")
                elif keys:
                    vals = l.strip("\n").split("\t")
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
    
    
    if len(self.picard_HsMetrics_data) > 0:
        
        # Write parsed data to a file
        self.write_data_file(self.picard_HsMetrics_data, 'multiqc_picard_HsMetrics')
        
        # Add to general stats table
        # Swap question marks with -1
        data = self.picard_HsMetrics_data
        for s_name in data:
            if data[s_name]['FOLD_ENRICHMENT'] == '?':
                data[s_name]['FOLD_ENRICHMENT'] = -1
                
        self.general_stats_headers['FOLD_ENRICHMENT'] = {
            'title': 'Fold Enrichment',
            'min': 0,
            'format': '{:.0f}',
            'scale': 'Blues',
        }
        self.general_stats_headers['PCT_TARGET_BASES_30X'] = {
            'title': 'Target Bases 30X',
            'description': 'Percent of target bases with coverage &ge; 30X',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'format': '{:.0f}%',
            'scale': 'RdYlGn',
            'modify': lambda x: self.multiply_hundred(x)
        }
        for s_name in data:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update( data[s_name] )

    table_html = plots.table.plot(data, _get_headers(data))
    if not isinstance(self.sections, list):
        self.sections = list()
    self.sections.append({
            'name': 'HSMetrics',
            'anchor': 'picard_hsmetrics',
            'content': table_html})
    # Return the number of detected samples to the parent module
    return len(self.picard_HsMetrics_data)


def _get_headers(data):
    header = {}
    for s in data:
        for h in data[s]:
            try:
                float(data[s][h])
            except:
                continue
            if h not in header:
                this = {
                'title': h.replace("_", " ")
                }
                if h.startswith("PCT"):
                    this.update({'max': 100, 'min': 0, 'format': '{:.0f}%', 'scale': 'RdYlGn', 'modify': lambda x: x * 100.0})
                header.update({h : this})
    return header
