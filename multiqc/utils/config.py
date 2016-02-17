#!/usr/bin/env python

""" MultiQC config module. Holds a single copy of
config variables to be used across all other modules """

from collections import defaultdict, OrderedDict
from datetime import datetime
import inspect
import logging
import os
import pkg_resources
import random
import sys
import yaml

import multiqc
from multiqc import (logger)
from multiqc.utils.log import init_log, LEVELS

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
data_dir_name = 'multiqc_report_data'
make_data_dir = True
data_format = 'tsv'
data_format_extensions = {'tsv': 'txt', 'json': 'json', 'yaml': 'yaml'}
fn_clean_exts = [ '.gz', '.fastq', '.fq', '.bam', '.sam', '.sra', '_tophat', '_star_aligned', '_fastqc' ]
fn_ignore_files = ['.DS_Store']
report_id = 'mqc_report_{}'.format(''.join(random.sample('abcdefghijklmnopqrstuvwxyz0123456789', 20)))
num_datasets_plot_limit = 50
log_filesize_limit = 1000000

#######################
# Module fn search patterns
#######################
# Load here so that they can be overwritten by user configs
searchp_fn = os.path.join( MULTIQC_DIR, 'utils', 'search_patterns.yaml')
with open(searchp_fn) as f:
    sp = yaml.load(f)

#######################
# Available modules
#######################
# Modules must be listed in setup.py under entry_points['multiqc.modules.v1']
# We load them here so that they can be overwritten by a config file if desired.

# Order that modules should appear in report. Try to list in order of analysis,
# eg. FastQC is usually the first step, so should be last in this list.
module_order = [
    # Post-alignment analysis results
    'qualimap', 'featureCounts', 'picard', 'preseq', 'samblaster',
    # Alignment tool stats
    'bismark', 'star', 'tophat', 'bowtie2', 'bowtie1',
    # Pre-alignment QC
    'cutadapt', 'fastq_screen', 'fastqc', 'skewer'
]

# Get all modules, including those from other extension packages
all_avail_modules = {}
avail_modules = OrderedDict()
for entry_point in pkg_resources.iter_entry_points('multiqc.modules.v1'):
    nicename = str(entry_point).split('=')[0].strip()
    all_avail_modules[nicename] = entry_point

# Start with modules not described above - probably plugins
for m in all_avail_modules.keys():
    if m not in module_order:
        avail_modules[m] = all_avail_modules[m]
        logger.debug("Module missing from order declaration: {}".format(m))

# Add known modules, in order defined above
for m in module_order:
    if m in all_avail_modules.keys():
        avail_modules[m] = all_avail_modules[m]


#######################
# Available templates
#######################
# Templates must be listed in setup.py under entry_points['multiqc.templates.v1']

# Get all templates, including those from other extension packages
avail_templates = {}
for entry_point in pkg_resources.iter_entry_points('multiqc.templates.v1'):
    nicename = str(entry_point).split('=')[0].strip()
    avail_templates[nicename] = entry_point

# Which template to use by default?
template = 'default'

#######################
# Check we have modules & templates
#######################
# Check that we were able to find some modules and templates
# If not, package probably hasn't been installed properly.
# Need to do this before click, else will throw cryptic error
# Note: Can't use logger here, not yet installed.
if len(avail_modules) == 0 or len(avail_templates) == 0:
    if len(avail_modules) == 0:
        print("Error - No MultiQC modules found.")
    if len(avail_templates) == 0:
        print("Error - No MultiQC templates found.")
    print("Could not load MultiQC - has it been installed? \n\
        Please either install with pip (pip install multiqc) or by using \n\
        the installation script (python setup.py install)")
    sys.exit(1)

#######################
# Overwrite defaults with config files
#######################
# Load and parse installation config file if we find it
try:
    yaml_config = os.path.join( os.path.dirname(MULTIQC_DIR), 'multiqc_config.yaml')
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

