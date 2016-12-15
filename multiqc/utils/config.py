#!/usr/bin/env python

""" MultiQC config module. Holds a single copy of
config variables to be used across all other modules """

from __future__ import print_function
from collections import OrderedDict
from datetime import datetime
import inspect
import os
import pkg_resources
import random
import subprocess
import sys
import yaml

import multiqc
from multiqc import logger

# Get the MultiQC version
version = pkg_resources.get_distribution("multiqc").version
cwd = os.getcwd()
script_path = os.path.dirname(os.path.realpath(__file__))
git_hash = None
git_hash_short = None
try:
    os.chdir(script_path)
    git_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'], stderr=subprocess.STDOUT)
    git_hash_short = git_hash[:7]
    version = '{} ({})'.format(version, git_hash_short)
except subprocess.CalledProcessError:
    pass
os.chdir(cwd)

# Constants
MULTIQC_DIR = os.path.dirname(os.path.realpath(inspect.getfile(multiqc)))

##### MultiQC Defaults
# Default MultiQC config
searchp_fn = os.path.join( MULTIQC_DIR, 'utils', 'config_defaults.yaml')
with open(searchp_fn) as f:
    configs = yaml.load(f)
    for c, v in configs.items():
        globals()[c] = v
# Module filename search patterns
searchp_fn = os.path.join( MULTIQC_DIR, 'utils', 'search_patterns.yaml')
with open(searchp_fn) as f:
    sp = yaml.load(f)

# Other defaults that can't be set in YAML
modules_dir = os.path.join(MULTIQC_DIR, 'modules')
creation_date = datetime.now().strftime("%Y-%m-%d, %H:%m")
working_dir = os.getcwd()
analysis_dir = [os.getcwd()]
output_dir = os.path.realpath(os.getcwd())
report_id = 'mqc_report_{}'.format(''.join(random.sample('abcdefghijklmnopqrstuvwxyz0123456789', 20)))

##### Available modules
# Modules must be listed in setup.py under entry_points['multiqc.modules.v1']
# Get all modules, including those from other extension packages
avail_modules = dict()
for entry_point in pkg_resources.iter_entry_points('multiqc.modules.v1'):
    nicename = str(entry_point).split('=')[0].strip()
    avail_modules[nicename] = entry_point

##### Available templates
# Templates must be listed in setup.py under entry_points['multiqc.templates.v1']
# Get all templates, including those from other extension packages
avail_templates = {}
for entry_point in pkg_resources.iter_entry_points('multiqc.templates.v1'):
    nicename = str(entry_point).split('=')[0].strip()
    avail_templates[nicename] = entry_point

##### Check we have modules & templates
# Check that we were able to find some modules and templates
# If not, package probably hasn't been installed properly.
# Need to do this before click, else will throw cryptic error
# Note: Can't use logger here, not yet initiated.
if len(avail_modules) == 0 or len(avail_templates) == 0:
    if len(avail_modules) == 0:
        print("Error - No MultiQC modules found.", file=sys.stderr)
    if len(avail_templates) == 0:
        print("Error - No MultiQC templates found.", file=sys.stderr)
    print("Could not load MultiQC - has it been installed? \n\
        Please either install with pip (pip install multiqc) or by using \n\
        the installation script (python setup.py install)", file=sys.stderr)
    sys.exit(1)

# Functions to load user config files. These are called by the main MultiQC script.
def mqc_load_userconfig(path=None):
    """ Overwrite config defaults with user config files """

    # Load and parse installation config file if we find it
    mqc_load_config(os.path.join( os.path.dirname(MULTIQC_DIR), 'multiqc_config.yaml'))

    # Load and parse a user config file if we find it
    mqc_load_config(os.path.expanduser('~/.multiqc_config.yaml'))

    # Load and parse a config file in this working directory if we find it
    mqc_load_config('multiqc_config.yaml')

    # Custom command line config
    if path is not None:
        mqc_load_config(path)


def mqc_load_config(yaml_config):
    """ Load and parse a config file if we find it """
    if os.path.isfile(yaml_config):
        try:
            with open(yaml_config) as f:
                new_config = yaml.load(f)
                logger.debug("Loading config settings from: {}".format(yaml_config))
                for c, v in new_config.items():
                    if c == 'sp':
                        # Merge filename patterns instead of replacing
                        sp.update(v)
                        logger.debug("Added to filename patterns: {}".format(sp))
                    elif c == 'extra_fn_clean_exts':
                        # Prepend to filename cleaning patterns instead of replacing
                        fn_clean_exts[0:0] = v
                        logger.debug("Added to filename clean extensions. Now looking for: {}".format(fn_clean_exts))
                    elif c == 'extra_fn_clean_trim':
                        # Prepend to filename cleaning patterns instead of replacing
                        fn_clean_trim[0:0] = v
                        logger.debug("Added to filename clean trimmings. Now looking for: {}".format(fn_clean_trim))
                    else:
                        logger.debug("New config '{}': {}".format(c, v))
                        globals()[c] = v
        except (IOError, AttributeError) as e:
            logger.debug("Config error: {}".format(e))
        except yaml.scanner.ScannerError as e:
            logger.error("Error parsing config YAML: {}".format(e))
            sys.exit(1)
    else:
        logger.debug("No MultiQC config found: {}".format(yaml_config))

