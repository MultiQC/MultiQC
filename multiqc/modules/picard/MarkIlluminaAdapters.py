#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard MarkIlluminaAdapters """

from collections import OrderedDict
import logging
import math
import os
import re

from multiqc import config

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find Picard MarkIlluminaAdapters reports and parse their data """

    # Set up vars
    self.picard_alignment_metrics = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files("picard/markilluminaadapters", filehandles=True):
        parsed_data = dict()
        s_name = None
        keys = None
