#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard RnaSeqMetrics """

from collections import OrderedDict
import logging
import os
import re

from multiqc.plots import linegraph, bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find Picard RnaSeqMetrics reports and parse their data """

    # Set up vars
    self.picard_RnaSeqMetrics_data = dict()
    self.picard_RnaSeqMetrics_histogram = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files('picard/rnaseqmetrics', filehandles=True):
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
            if 'rnaseqmetrics' in l.lower() and 'INPUT' in l:
                s_name = None
                # Pull sample name from input
                fn_search = re.search("INPUT=\[?([^\\s]+)\]?", l)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1))
                    s_name = self.clean_s_name(s_name, f['root'])

            if s_name is not None:
                if 'rnaseqmetrics' in l.lower() and '## METRICS CLASS' in l:
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
                    # Calculate some extra numbers
                    if 'PF_BASES' in keys and 'PF_ALIGNED_BASES' in keys:
                        self.picard_RnaSeqMetrics_data[s_name]['PF_NOT_ALIGNED_BASES'] = \
                            self.picard_RnaSeqMetrics_data[s_name]['PF_BASES'] - self.picard_RnaSeqMetrics_data[s_name]['PF_ALIGNED_BASES']

            if s_name is not None and 'normalized_position	All_Reads.normalized_coverage' in l:
                self.picard_RnaSeqMetrics_histogram[s_name] = dict()
                in_hist = True

        for key in list(self.picard_RnaSeqMetrics_data.keys()):
            if len(self.picard_RnaSeqMetrics_data[key]) == 0:
                self.picard_RnaSeqMetrics_data.pop(key, None)
        for s_name in list(self.picard_RnaSeqMetrics_histogram.keys()):
            if len(self.picard_RnaSeqMetrics_histogram[s_name]) == 0:
                self.picard_RnaSeqMetrics_histogram.pop(s_name, None)
                log.debug("Ignoring '{}' histogram as no data parsed".format(s_name))

    # Filter to strip out ignored sample names
    self.picard_RnaSeqMetrics_data = self.ignore_samples(self.picard_RnaSeqMetrics_data)

    if len(self.picard_RnaSeqMetrics_data) > 0:

        # Write parsed data to a file
        self.write_data_file(self.picard_RnaSeqMetrics_data, 'multiqc_picard_RnaSeqMetrics')

        # Add to general stats table
        GenStatsHeaders = OrderedDict()
        GenStatsHeaders['PCT_RIBOSOMAL_BASES'] = {
            'title': '% rRNA',
            'description': 'Percent of aligned bases overlapping ribosomal RNA regions',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Reds',
        }
        GenStatsHeaders['PCT_MRNA_BASES'] = {
            'title': '% mRNA',
            'description': 'Percent of aligned bases overlapping UTRs and coding regions of mRNA transcripts',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Greens',
        }
        self.general_stats_addcols(self.picard_RnaSeqMetrics_data, GenStatsHeaders)

        # Bar plot of bases assignment
        bg_cats = OrderedDict()
        bg_cats['CODING_BASES'] = { 'name': 'Coding' }
        bg_cats['UTR_BASES'] = { 'name': 'UTR' }
        bg_cats['INTRONIC_BASES'] = { 'name': 'Intronic' }
        bg_cats['INTERGENIC_BASES'] = { 'name': 'Intergenic' }
        bg_cats['RIBOSOMAL_BASES'] = { 'name': 'Ribosomal' }
        bg_cats['PF_NOT_ALIGNED_BASES'] = { 'name': 'PF not aligned' }
        self.add_section (
            name = 'RnaSeqMetrics Assignment',
            anchor = 'picard-rna-assignment',
            description = 'Number of bases in primary alignments that align to regions in the reference genome.',
            plot = bargraph.plot(self.picard_RnaSeqMetrics_data, bg_cats)
        )

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
            self.add_section (
                name = 'Gene Coverage',
                anchor = 'picard-rna-coverage',
                plot = linegraph.plot(self.picard_RnaSeqMetrics_histogram, pconfig)
            )


    # Return the number of detected samples to the parent module
    return len(self.picard_RnaSeqMetrics_data)

