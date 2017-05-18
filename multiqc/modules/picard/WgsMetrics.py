#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard WgsMetrics """

from collections import OrderedDict
import logging
import os
import re

from multiqc import config
from multiqc.plots import linegraph, bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find Picard WgsMetrics reports and parse their data """

    # Set up vars
    self.picard_wgsmetrics_data = dict()
    self.picard_wgsmetrics_histogram = dict()
    self.picard_wgsmetrics_samplestats = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files('picard/wgs_metrics', filehandles=True):
        s_name = None
        in_hist = False
        for l in f['f']:

            # Catch the histogram values
            if s_name is not None and in_hist is True:
                try:
                    sections = l.split("\t")
                    cov = int(sections[0])
                    count = int(sections[1])
                    self.picard_wgsmetrics_histogram[s_name][cov] = count
                except ValueError:
                    # Reset in case we have more in this log file
                    s_name = None
                    in_hist = False

            # New log starting
            if 'WgsMetrics' in l and 'INPUT' in l:
                s_name = None
                # Pull sample name from input
                fn_search = re.search("INPUT=\[?([^\\s]+)\]?", l)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1))
                    s_name = self.clean_s_name(s_name, f['root'])

            if s_name is not None:
                if 'CollectWgsMetrics$WgsMetrics' in l and '## METRICS CLASS' in l:
                    if s_name in self.picard_wgsmetrics_data:
                        log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
                    self.add_data_source(f, s_name, section='WgsMetrics')
                    self.picard_wgsmetrics_data[s_name] = dict()
                    keys = f['f'].readline().strip("\n").split("\t")
                    vals = f['f'].readline().strip("\n").split("\t")
                    if len(vals) == len(keys):
                        for i, k in enumerate(keys):
                            try:
                                self.picard_wgsmetrics_data[s_name][k] = float(vals[i])
                            except ValueError:
                                self.picard_wgsmetrics_data[s_name][k] = vals[i]

                    # Skip lines on to histogram
                    next(f['f'])
                    next(f['f'])
                    next(f['f'])

                    self.picard_wgsmetrics_histogram[s_name] = OrderedDict()
                    in_hist = True

        for key in list(self.picard_wgsmetrics_data.keys()):
            if len(self.picard_wgsmetrics_data[key]) == 0:
                self.picard_wgsmetrics_data.pop(key, None)
        for s_name in list(self.picard_wgsmetrics_histogram.keys()):
            if len(self.picard_wgsmetrics_histogram[s_name]) == 0:
                self.picard_wgsmetrics_histogram.pop(s_name, None)
                log.debug("Ignoring '{}' histogram as no data parsed".format(s_name))

    # Filter to strip out ignored sample names
    self.picard_wgsmetrics_data = self.ignore_samples(self.picard_wgsmetrics_data)

    if len(self.picard_wgsmetrics_data) > 0:

        # Write parsed data to a file
        self.write_data_file(self.picard_wgsmetrics_data, 'multiqc_picard_wgsmetrics')

        # Add to general stats table
        self.general_stats_headers['MEDIAN_COVERAGE'] = {
            'title': 'Median Coverage',
            'description': 'The median coverage in bases of the genome territory, after all filters are applied.',
            'min': 0,
            'suffix': 'X',
            'scale': 'GnBu',
        }

        # user configurable coverage level
        try:
            covs = config.picard_config['general_stats_target_coverage']
            assert type(covs) == list
            assert len(covs) > 0
            covs = [str(i) for i in covs]
            log.debug("Custom Picard coverage thresholds: {}".format(", ".join([i for i in covs])))
        except (AttributeError, TypeError, AssertionError):
            covs = ['30']
        for c in covs:
            self.general_stats_headers['PCT_{}X'.format(c)] = {
                'id': 'picard_target_bases_{}X'.format(c),
                'title': 'Bases &ge; {}X'.format(c),
                'description': 'Percent of target bases with coverage &ge; {}X'.format(c),
                'max': 100,
                'min': 0,
                'suffix': '%',
                'format': '{:,.0f}',
                'scale': 'RdYlGn',
                'modify': lambda x: self.multiply_hundred(x)
            }

        for s_name in self.picard_wgsmetrics_data:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update( self.picard_wgsmetrics_data[s_name] )

        # Section with histogram plot
        if len(self.picard_wgsmetrics_histogram) > 0:

            # Figure out where to cut histogram tail
            max_cov = 10
            for s_name, samp in self.picard_wgsmetrics_histogram.items():
                total = float( sum( samp.values() ) )
                running_total = 0
                for k, v in samp.items():
                    running_total += v
                    if running_total > total * 0.99:
                        max_cov = max(k, max_cov)
                        break

            # Cut histogram tail and make a normalised percentage version of the data plus dropoff
            data = {}
            data_percent = {}
            maxval = 0
            for s_name, samp in self.picard_wgsmetrics_histogram.items():
                data[s_name] = OrderedDict()
                data_percent[s_name] = OrderedDict()
                total = float( sum( samp.values() ) )
                cumulative = 0
                for k, v in samp.items():
                    if k <= max_cov:
                        cumulative += v
                        data[s_name][k] = v
                        maxval = max(maxval, v)
                        data_percent[s_name][k] = 100 - (cumulative/total)*100
                    else:
                        break

            # Plot the data and add section
            pconfig = {
                'id': 'picard_wgs_metrics_histogram',
                'title': 'Picard: WGS Coverage',
                'ylab': 'Percentage of Bases',
                'xlab': 'Fold Coverage',
                'xDecimals': False,
                'tt_label': '<b>{point.x} X</b>: {point.y:.1f}',
                'ymin': 0,
                'ymax': 100,
                'data_labels': [
                    {'name': 'Percentage Drop-Off', 'ylab': 'Percentage of Bases', 'ymax': 100},
                    {'name': 'Counts Histogram', 'ylab': 'Coverage', 'ymax': maxval}
                ]
            }
            self.add_section (
                name = 'WGS Coverage',
                anchor = 'picard-wgsmetrics-cov',
                description = 'The number of bases in the genome territory for each fold coverage. ' +
                            'Note that final 1% of data is hidden to prevent very long tails.',
                plot = linegraph.plot([data_percent, data], pconfig)
            )

            # Bar plot of ignored bases
            pdata = dict()
            for s_name, data in self.picard_wgsmetrics_data.items():
                pdata[s_name] = dict()
                pdata[s_name]['PCT_EXC_MAPQ'] = data['PCT_EXC_MAPQ'] * 100.0
                pdata[s_name]['PCT_EXC_DUPE'] = data['PCT_EXC_DUPE'] * 100.0
                pdata[s_name]['PCT_EXC_UNPAIRED'] = data['PCT_EXC_UNPAIRED'] * 100.0
                pdata[s_name]['PCT_EXC_BASEQ'] = data['PCT_EXC_BASEQ'] * 100.0
                pdata[s_name]['PCT_EXC_OVERLAP'] = data['PCT_EXC_OVERLAP'] * 100.0
                pdata[s_name]['PCT_EXC_CAPPED'] = data['PCT_EXC_CAPPED'] * 100.0

            keys = OrderedDict()
            keys['PCT_EXC_MAPQ'] = {'name': 'Low mapping quality'}
            keys['PCT_EXC_DUPE'] = {'name': 'Duplicates reads'}
            keys['PCT_EXC_UNPAIRED'] = {'name': 'No mapped mate pair'}
            keys['PCT_EXC_BASEQ'] = {'name': 'Low base quality'}
            keys['PCT_EXC_OVERLAP'] = {'name': 'Overlapping insert'}
            keys['PCT_EXC_CAPPED'] = {'name': 'Over capped coverage'}

            # Config for the plot
            pconfig = {
                'id': 'picard_wgs_metrics_bases',
                'title': 'Picard: WGS Filtered Bases',
                'cpswitch': False,
                'ylab': '% Bases',
                'ymax': 100
            }

            self.add_section (
                name = 'WGS Filtered Bases',
                anchor = 'picard-wgsmetrics-bases',
                description = 'For more information about the filtered categories, see the '+
                           '<a href="http://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectWgsMetrics.WgsMetrics" target="_blank">Picard documentation</a>.',
                plot = bargraph.plot(pdata, keys, pconfig)
            )


    # Return the number of detected samples to the parent module
    return len(self.picard_wgsmetrics_data)

