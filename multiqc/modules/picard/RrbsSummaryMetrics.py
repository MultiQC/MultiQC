#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard RrbsSummaryMetrics """

from collections import OrderedDict
import logging
import os
import re

from multiqc import config
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)

def parse_reports(self):
    """ Find Picard RrbsSummaryMetrics reports and parse their data """

    # Set up vars
    self.picard_rrbs_metrics = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files(config.sp['picard']['rrbs_metrics'], filehandles=True):
        parsed_data = dict()
        s_name = None
        keys = None
        for l in f['f']:
            # New log starting
            if 'CollectRrbsMetrics' in l and 'INPUT' in l:
                s_name = None
                keys = None
                # Pull sample name from input
                fn_search = re.search("INPUT=\[?([^\\s]+)\]?", l)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1))
                    s_name = self.clean_s_name(s_name, f['root'])
                    parsed_data[s_name] = dict()

            if s_name is not None:
                if 'RrbsSummaryMetrics' in l and '## METRICS CLASS' in l:
                    keys = f['f'].readline().strip("\n").split("\t")
                elif keys:
                    vals = l.strip("\n").split("\t")
                    if len(vals) == len(keys):
                        for i, k in enumerate(keys):
                            try:
                                parsed_data[s_name][k] = float(vals[i])
                            except ValueError:
                                parsed_data[s_name][k] = vals[i]
                    else:
                        s_name = None
                        keys = None

        # Remove empty dictionaries
        for s_name in list(parsed_data.keys()):
            if len(parsed_data[s_name]) == 0:
                parsed_data.pop(s_name, None)

        # Collect parsed data
        for s_name in parsed_data.keys():
            if s_name in self.picard_rrbs_metrics:
                log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
            self.add_data_source(f, s_name, section='RrbsSummaryMetrics')
            self.picard_rrbs_metrics[s_name] = parsed_data[s_name]

    if len(self.picard_rrbs_metrics) > 0:

        # Write parsed data to a file
        self.write_data_file(self.picard_rrbs_metrics, 'multiqc_picard_RrbsSummaryMetrics')

        # Add to general stats table
        self.general_stats_headers['PCT_CPG_BASES_CONVERTED'] = {
            'title': '% CpG Methylated',
            'description': 'Percentage of times a CpG cytosine was converted',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'format': '{:.0f}%',
            'scale': 'RdYlGn-rev',
            'modify': lambda x: 100 - self.multiply_hundred(x)
        }
        self.general_stats_headers['PCT_NON_CPG_BASES_CONVERTED'] = {
            'title': '% Non-CpG Methylated',
            'description': 'Percentage of times a non-CpG cytosine was converted',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'format': '{:.0f}%',
            'scale': 'RdYlGn',
            'modify': lambda x: 100 - self.multiply_hundred(x)
        }
        self.general_stats_headers['MEDIAN_CPG_COVERAGE'] = {
            'title': 'Median CpG Cov',
            'description': 'Median coverage of CpG sites',
            'min': 0
        }
        for s_name in self.picard_rrbs_metrics:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update( self.picard_rrbs_metrics[s_name] )



        # Make the bar plot of alignment read count
        pdata_cpg = dict()
        pdata_noncpg = dict()
        for s_name in self.picard_rrbs_metrics.keys():
            pdata_cpg[s_name] = dict()
            pdata_cpg[s_name]['converted'] = self.picard_rrbs_metrics[s_name]['CPG_BASES_CONVERTED']
            pdata_cpg[s_name]['not_converted'] = self.picard_rrbs_metrics[s_name]['CPG_BASES_SEEN'] - self.picard_rrbs_metrics[s_name]['CPG_BASES_CONVERTED']
            pdata_noncpg[s_name] = dict()
            pdata_noncpg[s_name]['converted'] = self.picard_rrbs_metrics[s_name]['NON_CPG_BASES']
            pdata_noncpg[s_name]['not_converted'] = self.picard_rrbs_metrics[s_name]['NON_CPG_BASES'] - self.picard_rrbs_metrics[s_name]['NON_CPG_CONVERTED_BASES']

        keys = OrderedDict()
        keys['not_converted'] = {'name': 'Unconverted Bases (Methylated)'}
        keys['converted'] = {'name': 'Converted Bases (Unmethylated)'}

        # Config for the plot
        pconfig = {
            'id': 'picard_rrbs_converted_bases_plot',
            'title': 'Picard: Converted Bases',
            'ylab': '# CpG Bases',
            'cpswitch_counts_label': 'Number of Bases',
            'data_labels': [
                {'name': 'CpG', 'ylab': '# CpG Bases'},
                {'name': 'Non-CpG', 'ylab': '# Non-CpG Bases'}
            ]
        }

        self.sections.append({
            'id': 'picard_rrbs_converted_bases',
            'name': 'RRBS Converted Bases',
            'anchor': 'picard-rrbssummary-convertedbases',
            'content': bargraph.plot([pdata_cpg, pdata_noncpg], [keys, keys], pconfig)
        })

    # Return the number of detected samples to the parent module
    return len(self.picard_rrbs_metrics)
