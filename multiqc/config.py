#!/usr/bin/env python

""" MultiQC config module. Holds a single copy of
config variables to be used across all other modules """

from collections import defaultdict, OrderedDict
from datetime import datetime
import inspect
import os

import multiqc

# Constants
VERSION = '0.1.1dev'
MULTIQC_DIR = os.path.dirname(os.path.realpath(inspect.getfile(multiqc)))

# Create the config variables, with defaults
title = None
prepend_dirs = False
creation_date = datetime.now().strftime("%Y-%m-%d, %H:%m")
working_dir = os.getcwd()
analysis_dir = os.getcwd()
output_dir = os.path.realpath(os.path.join(os.getcwd(), 'multiqc_report'))
template = 'default'
template_fn = os.path.join(output_dir, 'multiqc_report.html')
general_stats = {
    'headers': OrderedDict(),
    'rows': defaultdict(lambda:dict())
}

# Load and parse configuration file if we find it
try:
    yaml_config = os.path.join(MULTIQC_DIR, 'multiqc_config.yaml')
    with open(yaml_config) as f:
        config = yaml.load(f)
        for c, v in config.items():
            globals()[c] = v
except IOError:
    pass

# These config vars are imported by all modules and can be updated by anything.
# The main launcher (scripts/multiqc) overwrites some of these variables
# with what has been given to it on the command line.

