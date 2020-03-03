#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard MeanQualityByCycle"""

import logging
import os
import re

from multiqc.plots import linegraph
from .util import read_sample_name

# Initialise the logger
log = logging.getLogger(__name__)

def parse_reports(self):
    """ Find Picard QualityByCycleMetrics reports and parse their data """

    all_data = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files('picard/quality_by_cycle', filehandles=True):
        lines = iter(f['f'])

        # read through the header of the file to obtain the
        # sample name
        clean_fn = lambda n: self.clean_s_name(n, f['root'])
        s_name = read_sample_name(lines, clean_fn)
        if s_name is None:
            continue

        sample_data = dict()

        try:
            # skip to the histogram
            line = next(lines)
            while not line.startswith('## HISTOGRAM'):
                line = next(lines)

            # check the header
            line = next(lines)
            headers = line.strip().split("\t")
            if headers != ['CYCLE', 'MEAN_QUALITY']:
                continue

            # slurp the data
            line = next(lines).rstrip()
            while line:
                fields = line.split('\t')
                fields[0] = int(fields[0])
                fields[1] = float(fields[1])
                sample_data[fields[0]] = dict(zip(headers, fields))
                line = next(lines).rstrip()

        except StopIteration:
            pass

        # append the data
        if sample_data:
            all_data[s_name] = sample_data

    if not all_data:
        return 0

    # Write parsed data to a file
    self.write_data_file(all_data, 'multiqc_picard_quality_by_cycle')

    # Plot the data and add section
    pconfig = {
        'id': 'picard_quality_by_cycle',
        'title': 'Picard: Mean Base Quality by Cycle',
        'ylab': 'Mean Base Quality',
        'xlab': 'Cycle #',
        'xDecimals': False,
        'tt_label': '<b>cycle {point.x}</b>: {point.y:.2f} %',
        'ymax': 100,
        'ymin': 0,
    }

    lg = {}
    for s_name in all_data:
        lg[s_name] = dict((cycle, data['MEAN_QUALITY']) for cycle, data in all_data[s_name].items())

    self.add_section (
        name = 'Mean Base Quality by Cycle',
        anchor = 'picard-quality-by-cycle',
        description = 'Plot shows the mean base quality by cycle.',
        plot = linegraph.plot([lg], pconfig)
    )

    # Return the number of detected samples to the parent module
    return len(all_data)
