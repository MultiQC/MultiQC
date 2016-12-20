#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard RnaSeqMetrics """

from collections import OrderedDict
import logging
import os
import re

from multiqc import config, plots

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find Picard RnaSeqMetrics reports and parse their data """

    # Set up vars
    self.picard_RnaSeqMetrics_data = dict()
    self.picard_RnaSeqMetrics_histogram = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files(config.sp['picard']['rnaseqmetrics'], filehandles=True):
        s_name = None
        in_hist = False
        for l in f['f']:
            # Catch the histogram values
            if s_name is not None and in_hist is True:
                try:
                    sections = l.split("\t")
                    pos = int(sections[0])
                    coverage = float(sections[1])
                    self.picard_RnaSeqMetrics_histogram[s_name][pos] = coverage
                except ValueError:
                    # Reset in case we have more in this log file
                    s_name = None
                    in_hist = False

            # New log starting
            if 'RnaSeqMetrics' in l and 'INPUT' in l:
                s_name = None
                # Pull sample name from input
                fn_search = re.search("INPUT=\[?([^\\s]+)\]?", l)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1))
                    s_name = self.clean_s_name(s_name, f['root'])

            if s_name is not None:
                if 'RnaSeqMetrics' in l and '## METRICS CLASS' in l:
                    if s_name in self.picard_RnaSeqMetrics_data:
                        log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
                    self.picard_RnaSeqMetrics_data[s_name] = dict()
                    self.picard_RnaSeqMetrics_histogram[s_name] = dict()
                    self.add_data_source(f, s_name, section='RnaSeqMetrics')
                    keys = f['f'].readline().strip("\n").split("\t")
                    vals = f['f'].readline().strip("\n").split("\t")
                    for i, k in enumerate(keys):
                        # Multiply percentages by 100
                        if k.startswith('PCT_'):
                            try:
                                vals[i] = float(vals[i]) * 100.0
                            except (ValueError, IndexError):
                                pass
                        # Save the key:value pairs
                        try:
                            self.picard_RnaSeqMetrics_data[s_name][k] = float(vals[i])
                        except ValueError:
                            self.picard_RnaSeqMetrics_data[s_name][k] = vals[i]
                        except IndexError:
                            pass # missing data

                    # Skip lines on to histogram
                    l = f['f'].readline()
                    l = f['f'].readline()
                    l = f['f'].readline().strip("\n")

                    self.picard_RnaSeqMetrics_histogram[s_name] = dict()
                    in_hist = True

        for key in list(self.picard_RnaSeqMetrics_data.keys()):
            if len(self.picard_RnaSeqMetrics_data[key]) == 0:
                self.picard_RnaSeqMetrics_data.pop(key, None)
        for s_name in list(self.picard_RnaSeqMetrics_histogram.keys()):
            if len(self.picard_RnaSeqMetrics_histogram[s_name]) == 0:
                self.picard_RnaSeqMetrics_histogram.pop(s_name, None)
                log.debug("Ignoring '{}' histogram as no data parsed".format(s_name))

    if len(self.picard_RnaSeqMetrics_data) > 0:

        # Write parsed data to a file
        self.write_data_file(self.picard_RnaSeqMetrics_data, 'multiqc_picard_RnaSeqMetrics')

        # Add to general stats table
        self.general_stats_headers['PCT_CODING_BASES'] = {
            'title': '% Coding',
            'description': 'Percent of coding bases',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'format': '{:.1f}%',
            'scale': 'Greens',
        }
        self.general_stats_headers['PCT_UTR_BASES'] = {
            'title': '% UTR',
            'description': 'Percent of UTR bases',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'format': '{:.1f}%',
            'scale': 'GnBu',
        }
        self.general_stats_headers['PCT_INTRONIC_BASES'] = {
            'title': '% Intronic',
            'description': 'Percent of intronic bases',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'format': '{:.1f}%',
            'scale': 'Reds',
        }
        self.general_stats_headers['PCT_INTERGENIC_BASES'] = {
            'title': '% Intergenic',
            'description': 'Percent of intergenic bases',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'format': '{:.1f}%',
            'scale': 'Reds',
        }
        self.general_stats_headers['MEDIAN_CV_COVERAGE'] = {
            'title': 'CV Cov',
            'description': 'Median CV coverage',
            'min': 0,
            'format': '{:.1f}',
            'scale': 'GnBu',
        }
        self.general_stats_headers['MEDIAN_5PRIME_BIAS'] = {
            'title': "5' bias",
            'description': "Median 5' bias",
            'min': 0,
            'format': '{:.2f}',
            'scale': 'GnBu',
        }
        self.general_stats_headers['MEDIAN_3PRIME_BIAS'] = {
            'title': "3' bias",
            'description': "Median 3' bias",
            'min': 0,
            'format': '{:.2f}',
            'scale': 'GnBu',
        }
        for s_name in self.picard_RnaSeqMetrics_data:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update( self.picard_RnaSeqMetrics_data[s_name] )

        # Section with histogram plot
        if len(self.picard_RnaSeqMetrics_histogram) > 0:
            # Plot the data and add section
            pconfig = {
                'smooth_points': 500,
                'smooth_points_sumcounts': [True, False],
                'id': 'picard_rna_coverage',
                'title': 'Normalized Coverage',
                'ylab': 'Coverage',
                'xlab': 'Percent through gene',
                'xDecimals': False,
                'tt_label': '<b>{point.x}%</b>: {point.y:.0f}',
                'ymin': 0,
            }
            self.sections.append({
                'id': 'picard_gene_coverage',
                'name': 'Gene Coverage',
                'anchor': 'picard-rna-coverage',
                'content': plots.linegraph.plot(self.picard_RnaSeqMetrics_histogram, pconfig)
            })


    # Return the number of detected samples to the parent module
    return len(self.picard_RnaSeqMetrics_data)

