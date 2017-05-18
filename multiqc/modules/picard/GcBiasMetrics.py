#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard InsertSizeMetrics """

import logging
import os
import re

from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find Picard InsertSizeMetrics reports and parse their data """

    # Set up vars
    self.picard_GCbias_data = dict()
    self.picard_GCbiasSummary_data = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files('picard/gcbias', filehandles=True):
        s_name = None
        gc_col = None
        cov_col = None
        for l in f['f']:
            # New log starting
            if 'GcBiasMetrics' in l and 'INPUT' in l:
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

                if 'GcBiasDetailMetrics' in l and '## METRICS CLASS' in l:
                    if s_name in self.picard_GCbias_data:
                        log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
                    self.add_data_source(f, s_name, section='GcBiasDetailMetrics')
                    self.picard_GCbias_data[s_name] = dict()
                    # Get header - find columns with the data we want
                    l = f['f'].readline()
                    s = l.strip("\n").split("\t")
                    gc_col = s.index('GC')
                    cov_col = s.index('NORMALIZED_COVERAGE')

                if 'GcBiasSummaryMetrics' in l and '## METRICS CLASS' in l:
                    if s_name in self.picard_GCbias_data:
                        log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
                    self.add_data_source(f, s_name, section='GcBiasSummaryMetrics')
                    self.picard_GCbiasSummary_data[s_name] = dict()

                    keys = f['f'].readline().rstrip("\n").split("\t")
                    vals = f['f'].readline().rstrip("\n").split("\t")
                    for i, k in enumerate(keys):
                        try:
                            self.picard_GCbiasSummary_data[s_name][k] = float(vals[i])
                        except ValueError:
                            self.picard_GCbiasSummary_data[s_name][k] = vals[i]


        for s_name in list(self.picard_GCbias_data.keys()):
            if len(self.picard_GCbias_data[s_name]) == 0:
                self.picard_GCbias_data.pop(s_name, None)
                log.debug("Removing {} as no data parsed".format(s_name))

        for s_name in list(self.picard_GCbiasSummary_data.keys()):
            if len(self.picard_GCbiasSummary_data[s_name]) == 0:
                self.picard_GCbiasSummary_data.pop(s_name, None)
                log.debug("Removing {} as no data parsed".format(s_name))


    # Filter to strip out ignored sample names
    self.picard_GCbias_data = self.ignore_samples(self.picard_GCbias_data)

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
        self.add_section (
            name = 'GC Coverage Bias',
            anchor = 'picard-gcbias',
            description = 'This plot shows bias in coverage across regions of the genome with varying GC content.'\
                ' A perfect library would be a flat line at <code>y = 1</code>.',
            plot = linegraph.plot(self.picard_GCbias_data, pconfig)
        )

    if len(self.picard_GCbiasSummary_data) > 0:
        # Write parsed summary data to a file
        self.write_data_file(self.picard_GCbiasSummary_data, 'multiqc_picard_gcbias')


    # Return the number of detected samples to the parent module
    return len(self.picard_GCbias_data)
