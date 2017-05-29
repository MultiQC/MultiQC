#!/usr/bin/env python

""" MultiQC submodule to parse output from RSeQC inner_distance.py
http://rseqc.sourceforge.net/#inner-distance-py """

from collections import OrderedDict
import logging

from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find RSeQC inner_distance frequency reports and parse their data """

    # Set up vars
    self.inner_distance = dict()
    self.inner_distance_pct = dict()

    # Go through files and parse data
    for f in self.find_log_files('rseqc/inner_distance'):
        if f['s_name'] in self.inner_distance:
            log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
        self.add_data_source(f, section='inner_distance')
        #saving to temporary variable fro SE checking later
        parsed_data = OrderedDict()
        for l in f['f'].splitlines():
            s = l.split()
            try:
                avg_pos = (float(s[0]) +  float(s[1])) / 2.0
                parsed_data[avg_pos] = float(s[2])
            except:
                # Don't bother running through whole file if wrong
                break

        # Only add if we actually found something i,e it was PE data
        if len(parsed_data) > 0:
            self.inner_distance[f['s_name']] = parsed_data

    # Filter to strip out ignored sample names
    self.inner_distance = self.ignore_samples(self.inner_distance)

    if len(self.inner_distance) > 0:

        # Make a normalised percentage version of the data
        for s_name in self.inner_distance:
            self.inner_distance_pct[s_name] = OrderedDict()
            total = sum( self.inner_distance[s_name].values() )
            for k, v in self.inner_distance[s_name].items():
                self.inner_distance_pct[s_name][k] = (v/total)*100

        # Add line graph to section
        pconfig = {
            'id': 'rseqc_inner_distance_plot',
            'title': 'RSeQC: Inner Distance',
            'ylab': 'Counts',
            'xlab': "Inner Distance (bp)",
            'tt_label': "<strong>{point.x} bp</strong>: {point.y:.2f}",
            'data_labels': [
                {'name': 'Counts', 'ylab': 'Counts'},
                {'name': 'Percentages', 'ylab': 'Percentage'}
            ]
        }
        self.add_section (
            name = 'Inner Distance',
            anchor = 'rseqc-inner_distance',
            description = '<a href="http://rseqc.sourceforge.net/#inner-distance-py" target="_blank">Inner Distance</a>' \
                " calculates the inner distance" \
                " (or insert size) between two paired RNA reads." \
                " Note that this can be negative if fragments overlap.",
            plot = linegraph.plot([self.inner_distance, self.inner_distance_pct], pconfig)
        )

    # Return number of samples found
    return len(self.inner_distance)
