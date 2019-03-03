#!/usr/bin/env python

""" MultiQC config module. Holds a single copy of
config variables to be used across all other modules """

from __future__ import print_function
from datetime import datetime
import inspect
import collections
import os
import pkg_resources
import subprocess
import sys
import yaml

import multiqc

# Default logger will be replaced by caller
import logging
logger = logging.getLogger(__name__)

# Get the MultiQC version
version = pkg_resources.get_distribution("multiqc").version
short_version = pkg_resources.get_distribution("multiqc").version
script_path = os.path.dirname(os.path.realpath(__file__))
git_hash = None
git_hash_short = None
try:
    git_hash = subprocess.check_output( ['git', 'rev-parse', 'HEAD'],
                                        cwd=script_path,
                                        stderr=subprocess.STDOUT,
                                        universal_newlines=True ).strip()
    git_hash_short = git_hash[:7]
    version = '{} ({})'.format(version, git_hash_short)
except:
    pass

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
data_tmp_dir = '/tmp' # will be overwritten by core script
modules_dir = os.path.join(MULTIQC_DIR, 'modules')
creation_date = datetime.now().strftime("%Y-%m-%d, %H:%M")
working_dir = os.getcwd()
analysis_dir = [os.getcwd()]
output_dir = os.path.realpath(os.getcwd())
megaqc_access_token = os.environ.get('MEGAQC_ACCESS_TOKEN')

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

##### Functions to load user config files. These are called by the main MultiQC script.
# Note that config files are loaded in a specific order and values can overwrite each other.
def mqc_load_userconfig(paths=()):
    """ Overwrite config defaults with user config files """

    # Load and parse installation config file if we find it
    mqc_load_config(os.path.join( os.path.dirname(MULTIQC_DIR), 'multiqc_config.yaml'))

    # Load and parse a user config file if we find it
    mqc_load_config(os.path.expanduser('~/.multiqc_config.yaml'))

    # Load and parse a config file path set in an ENV variable if we find it
    if os.environ.get('MULTIQC_CONFIG_PATH') is not None:
        mqc_load_config( os.environ.get('MULTIQC_CONFIG_PATH') )

    # Load and parse a config file in this working directory if we find it
    mqc_load_config('multiqc_config.yaml')

    # Custom command line config
    for p in paths:
        mqc_load_config(p)


def mqc_load_config(yaml_config):
    """ Load and parse a config file if we find it """
    if os.path.isfile(yaml_config):
        try:
            with open(yaml_config) as f:
                new_config = yaml.load(f)
                logger.debug("Loading config settings from: {}".format(yaml_config))
                mqc_add_config(new_config, yaml_config)
        except (IOError, AttributeError) as e:
            logger.debug("Config error: {}".format(e))
        except yaml.scanner.ScannerError as e:
            logger.error("Error parsing config YAML: {}".format(e))
            sys.exit(1)
    else:
        logger.debug("No MultiQC config found: {}".format(yaml_config))

def mqc_cl_config(cl_config):
    for clc_str in cl_config:
        try:
            parsed_clc = yaml.load(clc_str)
            # something:var fails as it needs a space. Fix this (a common mistake)
            if isinstance(parsed_clc, str) and ':' in clc_str:
                clc_str = ': '.join(clc_str.split(':'))
                parsed_clc = yaml.load(clc_str)
            assert(isinstance(parsed_clc, dict))
        except yaml.scanner.ScannerError as e:
            logger.error("Could not parse command line config: {}\n{}".format(clc_str, e))
        except AssertionError:
            logger.error("Could not parse command line config: {}".format(clc_str))
        else:
            logger.debug("Found command line config: {}".format(parsed_clc))
            mqc_add_config(parsed_clc)

def mqc_add_config(conf, conf_path=None):
    """ Add to the global config with given MultiQC config dict """
    global fn_clean_exts, fn_clean_trim
    for c, v in conf.items():
        if c == 'sp':
            # Merge filename patterns instead of replacing
            sp.update(v)
            logger.debug("Added to filename patterns: {}".format(v))
        elif c == 'extra_fn_clean_exts':
            # Prepend to filename cleaning patterns instead of replacing
            fn_clean_exts[0:0] = v
            logger.debug("Added to filename clean extensions: {}".format(v))
        elif c == 'extra_fn_clean_trim':
            # Prepend to filename cleaning patterns instead of replacing
            fn_clean_trim[0:0] = v
            logger.debug("Added to filename clean trimmings: {}".format(v))
        elif c in ['custom_logo'] and v:
            # Resolve file paths - absolute or cwd, or relative to config file
            fpath = v
            if os.path.exists(v):
                fpath = os.path.abspath(v)
            elif conf_path is not None and os.path.exists(os.path.join(os.path.dirname(conf_path), v)):
                fpath = os.path.abspath(os.path.join(os.path.dirname(conf_path), v))
            else:
                logger.error("Config '{}' path not found, skipping ({})".format(c, fpath))
                continue
            logger.debug("New config '{}': {}".format(c, fpath))
            update({c: fpath})
        else:
            logger.debug("New config '{}': {}".format(c, v))
            update({c: v})

#### Function to load file containing a list of alternative sample-name swaps
# Essentially a fancy way of loading stuff into the sample_names_rename config var
# As such, can also be done directly using a config file
def load_sample_names(snames_file):
    global sample_names_rename_buttons, sample_names_rename
    num_cols = None
    try:
        with open(snames_file) as f:
            logger.debug("Loading sample renaming config settings from: {}".format(snames_file))
            for l in f:
                s = l.strip().split("\t")
                if len(s) > 1:
                    # Check that we have consistent numbers of columns
                    if num_cols is None:
                        num_cols = len(s)
                    elif num_cols != len(s):
                        logger.warn("Inconsistent number of columns found in sample names file (skipping line): '{}'".format(l.strip()))
                    # Parse the line
                    if len(sample_names_rename_buttons) == 0:
                        sample_names_rename_buttons = s
                    else:
                        sample_names_rename.append(s)
                elif len(l.strip()) > 0:
                    logger.warn("Sample names file line did not have columns (must use tabs): {}".format(l.strip()))
    except (IOError, AttributeError) as e:
        logger.error("Error loading sample names file: {}".format(e))
    logger.debug("Found {} sample renaming patterns".format(len(sample_names_rename_buttons)))

def update(u):
    return update_dict(globals(), u)

def update_dict(d, u):
    """ Recursively updates nested dict d from nested dict u
    """
    for key, val in u.items():
        if isinstance(val, collections.Mapping):
            d[key] = update_dict(d.get(key, {}), val)
        else:
            d[key] = u[key]
    return d
