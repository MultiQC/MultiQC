#!/usr/bin/env python

""" MultiQC submodule to parse output from RSeQC read_duplication.py
http://rseqc.sourceforge.net/#read-duplication-py """

from collections import OrderedDict
import logging

from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find RSeQC read_duplication reports and parse their data """

    # Set up vars
    self.read_dups = dict()

    # Go through files and parse data
    for f in self.find_log_files('rseqc/read_duplication_pos'):
        if f['f'].startswith('Occurrence	UniqReadNumber'):
            if f['s_name'] in self.read_dups:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            self.add_data_source(f, section='read_duplication')
            self.read_dups[f['s_name']] = OrderedDict()
            for l in f['f'].splitlines():
                s = l.split()
                try:
                    if int(s[0]) <= 500:
                        self.read_dups[f['s_name']][int(s[0])] = int(s[1])
                except:
                    pass

    # Filter to strip out ignored sample names
    self.read_dups = self.ignore_samples(self.read_dups)

    if len(self.read_dups) > 0:

        # Add line graph to section
        pconfig = {
            'smooth_points': 200,
            'id': 'rseqc_read_dups_plot',
            'title': 'RSeQC: Read Duplication',
            'ylab': 'Number of Reads (log10)',
            'xlab': "Occurance of read",
            'yLog': True,
            'tt_label': "<strong>{point.x} occurances</strong>: {point.y} reads",
        }

        self.add_section (
            name = 'Read Duplication',
            anchor = 'rseqc-read_dups',
            description = '<a href="http://rseqc.sourceforge.net/#read-duplication-py" target="_blank">read_duplication.py</a>' \
                " calculates how many alignment positions have a certain number of exact duplicates."\
                " Note - plot truncated at 500 occurances and binned.</p>",
            plot = linegraph.plot(self.read_dups, pconfig)
        )

    # Return number of samples found
    return len(self.read_dups)



