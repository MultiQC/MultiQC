#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard MeanQualityByCycle"""

import logging
import os
import re

from multiqc.plots import linegraph
from .util import read_histogram

# Initialise the logger
log = logging.getLogger(__name__)

def parse_reports(self):
    """ Find Picard QualityByCycleMetrics reports and parse their data """

    headers = ['CYCLE', 'MEAN_QUALITY']
    formats = [int, float]
    all_data = read_histogram(self, 'picard/quality_by_cycle', 'MeanQualityByCycle', headers, formats)

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
