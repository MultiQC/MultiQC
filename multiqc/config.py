#!/usr/bin/env python

""" MultiQC config module. Holds a single copy of
config variables to be used across all other modules """

from collections import defaultdict, OrderedDict
from datetime import datetime
import inspect
import logging
import os
import pkg_resources
import yaml

import multiqc
from multiqc import (logger)
from multiqc.log import init_log, LEVELS

# Constants
MULTIQC_DIR = os.path.dirname(os.path.realpath(inspect.getfile(multiqc)))

#######################
# Main Config Variables
#######################
title = None
prepend_dirs = False
modules_dir = os.path.join(MULTIQC_DIR, 'modules')
creation_date = datetime.now().strftime("%Y-%m-%d, %H:%m")
working_dir = os.getcwd()
analysis_dir = [os.getcwd()]
output_dir = os.path.realpath(os.getcwd())
output_fn_name = 'multiqc_report.html'
data_dir_name = 'multiqc_data'
general_stats = {
    'headers': OrderedDict(),
    'rows': defaultdict(lambda:dict())
}
fn_clean_exts = [ '.gz', '.fastq', '.fq', '.bam', '.sam', '_tophat', '_star_aligned', '_fastqc' ]
fn_ignore_files = ['.DS_Store', 'glyphicons-halflings-regular.woff2']


#######################
# Available modules
#######################
# Modules must be listed in setup.py under entry_points['multiqc.modules.v1']

# Order that modules should appear in report. Try to list in order of analysis,
# eg. FastQC is usually the first step, so should be last in this list.
module_order = [
    # Post-alignment analysis results
    'qualimap', 'featureCounts', 'picard',
    # Alignment tool stats
    'bismark', 'star', 'tophat', 'bowtie2', 'bowtie1',
    # Pre-alignment QC
    'cutadapt', 'fastq_screen', 'fastqc'
]

# Get all modules, including those from other extension packages
all_avail_modules = {}
for entry_point in pkg_resources.iter_entry_points('multiqc.modules.v1'):
    nicename = str(entry_point).split('=')[0].strip()
    all_avail_modules[nicename] = entry_point

# Create ordered list of modules, as defined above
avail_modules = OrderedDict()
for m in module_order:
    if m in all_avail_modules.keys():
        avail_modules[m] = all_avail_modules[m]

# Add on any not described in the order above
for m in all_avail_modules.keys():
    if m not in module_order:
        avail_modules[m] = all_avail_modules[m]
        logger.debug("Module missing from order declaration: {}".format(m))


#######################
# Available templates
#######################
# Templates must be listed in setup.py under entry_points['multiqc.templates.v1']

# Get all templates, including those from other extension packages
template_entry_points = pkg_resources.iter_entry_points('multiqc.templates.v1')

avail_templates = [ d for d in os.listdir(os.path.join(MULTIQC_DIR, 'templates'))
                if os.path.isdir(os.path.join(MULTIQC_DIR, 'templates', d)) ]

# Which template to use by default?
template = 'default'
template_dir = os.path.join(MULTIQC_DIR, 'templates', template.strip())
template_fn = os.path.join(template_dir, 'multiqc_report.html')

#######################
# Overwrite defaults with config files
#######################
# Load and parse installation config file if we find it
try:
    yaml_config = os.path.join(MULTIQC_DIR, 'multiqc_config.yaml')
    with open(yaml_config) as f:
        config = yaml.load(f)
        for c, v in config.items():
            globals()[c] = v
except (IOError, AttributeError):
    pass

# Load and parse a user config file if we find it
try:
    yaml_config = os.path.expanduser('~/.multiqc_config.yaml')
    with open(yaml_config) as f:
        config = yaml.load(f)
        for c, v in config.items():
            globals()[c] = v
except (IOError, AttributeError):
    pass

# These config vars are imported by all modules and can be updated by anything.
# The main launcher (scripts/multiqc) overwrites some of these variables
# with what has been given to it on the command line.

