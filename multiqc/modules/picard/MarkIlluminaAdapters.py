#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard MarkIlluminaAdapters """

from collections import OrderedDict
import logging
import math
import os
import re

from multiqc import config
from multiqc.plots import linegraph
from .util import read_histogram

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find Picard MarkIlluminaAdapters reports and parse their data """

    headers = ["clipped_bases", "read_count"]
    formats = [int, int]
    all_data = read_histogram(
        self, "picard/markilluminaadapters", "MarkIlluminaAdapters", headers, formats
    )

    if not all_data:
        return 0

    # Write parsed data to a file
    self.write_data_file(all_data, "multiqc_picard_mark_illumina_adapters")

    # Plot the data and add section
    pconfig = {
        "id": "picard_mark_illumina_adapters",
        "title": "Picard: Mark Illumina Adapters",
        "ylab": "Clipped Bases",
        "xlab": "Read Number",
        "xDecimals": False,
        "tt_label": "<b>cycle {point.x}</b>: {point.y:.2f}",
        "ymin": 0,
    }

    lg = {}
    for s_name in all_data:
        lg[s_name] = {clipped_bases:data['read_count'] for clipped_bases, data in all_data[s_name].items()}

    self.add_section(
        name="Mark Illumina Adapters",
        description="Number of Clipped Bases by Read",
        anchor="picard-mark-illumina-adapters",
        plot=linegraph.plot(lg, pconfig),
    )

    # Return the number of detected samples to the parent module
    return len(all_data)
