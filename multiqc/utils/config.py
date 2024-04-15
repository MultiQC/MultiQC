""" MultiQC config module. Holds a single copy of
config variables to be used across all other modules """

from pathlib import Path
from typing import List, Dict, Optional, Union

import inspect

# Default logger will be replaced by caller
import logging
import os
import subprocess
import sys
from datetime import datetime

import importlib_metadata
import yaml
import pyaml_env

import multiqc
from multiqc.utils.util_functions import strtobool

logger = logging.getLogger("multiqc")

# Get the MultiQC version
version = importlib_metadata.version("multiqc")
short_version = version
git_hash = None
git_hash_short = None
script_path = str(Path(__file__).parent)  # dynamically used by util_functions.multiqc_dump_json()
git_root = None
try:
    git_root = subprocess.check_output(
        ["git", "rev-parse", "--show-toplevel"], cwd=script_path, stderr=subprocess.STDOUT, universal_newlines=True
    ).strip()
    git_root = Path(git_root)
    # .git
    # multiqc/
    #   utils/
    #       config.py  <- __file__
    expected_git_root = Path(script_path).parent.parent
    if git_root == expected_git_root:
        git_hash = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=script_path, stderr=subprocess.STDOUT, universal_newlines=True
        ).strip()
        git_hash_short = git_hash[:7]
        version = f"{version} ({git_hash_short})"
except Exception:
    pass

# Constants
MULTIQC_DIR = os.path.dirname(os.path.realpath(inspect.getfile(multiqc)))

##### MultiQC Defaults
# Declaring variables for static typing
title: str
subtitle: str
intro_text: str
report_comment: str
report_header_info: List[Dict[str, str]]
show_analysis_paths: bool
show_analysis_time: bool
config_file: str
custom_logo: str
custom_logo_url: str
custom_logo_title: str
custom_css_files: List[str]
simple_output: bool
template: str
profile_runtime: bool
pandoc_template: str
read_count_multiplier: float
read_count_prefix: str
read_count_desc: str
long_read_count_multiplier: float
long_read_count_prefix: str
long_read_count_desc: str
base_count_multiplier: float
base_count_prefix: str
base_count_desc: str
output_fn_name: str
data_dir_name: str
plots_dir_name: str
data_format: str
force: bool
no_ansi: bool
quiet: bool
prepend_dirs: bool
prepend_dirs_depth: int
prepend_dirs_sep: str
file_list: bool
require_logs: bool
version_check_url: str

make_data_dir: bool
zip_data_dir: bool
data_dump_file: bool
megaqc_url: str
megaqc_access_token: str
megaqc_timeout: float
export_plots: bool
make_report: bool
make_pdf: bool
plots_force_flat: bool
plots_force_interactive: bool
plots_flat_numseries: int
num_datasets_plot_limit: int
lineplot_style: str
lineplot_max_samples: int
barplot_legend_on_bottom: bool
violin_downsample_after: int
violin_min_threshold_outliers: int
violin_min_threshold_no_points: int
collapse_tables: bool
max_table_rows: int
table_columns_visible: Dict
table_columns_placement: Dict
table_columns_name: Dict
table_cond_formatting_colours: List[Dict[str, str]]
table_cond_formatting_rules: Dict[str, List[Dict[str, str]]]
decimalPoint_format: str
thousandsSep_format: str
remove_sections: List
section_comments: Dict
lint: bool  # Deprecated since v1.17
strict: bool
development: bool
custom_plot_config: Dict
custom_table_header_config: Dict
software_versions: Dict
ignore_symlinks: bool
ignore_images: bool
fn_ignore_dirs: List[str]
fn_ignore_paths: List[str]
sample_names_ignore: List[str]
sample_names_ignore_re: List[str]
sample_names_rename_buttons: List[str]
sample_names_replace: Dict
sample_names_replace_regex: bool
sample_names_replace_exact: bool
sample_names_replace_complete: bool
sample_names_rename: List
show_hide_buttons: List
show_hide_patterns: List
show_hide_regex: List
show_hide_mode: List
no_version_check: bool
log_filesize_limit: int
filesearch_lines_limit: int
report_readerrors: int
skip_generalstats: int
skip_versions_section: int
disable_version_detection: int
versions_table_group_header: str
data_format_extensions: Dict[str, str]
export_plot_formats: List[str]
filesearch_file_shared: List[str]
custom_content: Dict
fn_clean_sample_names: bool
use_filename_as_sample_name: bool
fn_clean_exts: List
fn_clean_trim: List
fn_ignore_files: List
top_modules: List
module_order: List[Union[str, Dict]]

# Populating the variables above from the default MultiQC config
config_defaults_path = os.path.join(MULTIQC_DIR, "utils", "config_defaults.yaml")
with open(config_defaults_path) as f:
    _default_config = yaml.safe_load(f)
for c, v in _default_config.items():
    globals()[c] = v

# Module filename search patterns
searchp_fn = os.path.join(MULTIQC_DIR, "utils", "search_patterns.yaml")
with open(searchp_fn) as f:
    sp = yaml.safe_load(f)

# Other defaults that can't be set in YAML
data_tmp_dir = "/tmp"  # will be overwritten by core script
modules_dir = os.path.join(MULTIQC_DIR, "modules")
creation_date = datetime.now().astimezone().strftime("%Y-%m-%d, %H:%M %Z")
working_dir = os.getcwd()
analysis_dir = [os.getcwd()]
output_dir = os.path.realpath(os.getcwd())
megaqc_access_token = os.environ.get("MEGAQC_ACCESS_TOKEN")

# Other variables that set only through the CLI
run_modules: List[str] = []
exclude_modules: List[str] = []
data_dir: Optional[str] = None
plots_tmp_dir: Optional[str] = None
plots_dir: Optional[str] = None
custom_data: Dict = {}
report_section_order: Dict = {}
output_fn: Optional[str] = None
megaqc_upload: bool = False

##### Available modules
# Modules must be listed in setup.py under entry_points['multiqc.modules.v1']
# Get all modules, including those from other extension packages
avail_modules = dict()
for entry_point in importlib_metadata.entry_points(group="multiqc.modules.v1"):
    nicename = entry_point.name
    avail_modules[nicename] = entry_point

##### Available templates
# Templates must be listed in setup.py under entry_points['multiqc.templates.v1']
# Get all templates, including those from other extension packages
avail_templates = {}
for entry_point in importlib_metadata.entry_points(group="multiqc.templates.v1"):
    nicename = entry_point.name
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
    print(
        "Could not load MultiQC - has it been installed? \n\
        Please either install with pip (pip install multiqc) or by using \n\
        the local files (pip install .)",
        file=sys.stderr,
    )
    sys.exit(1)


##### Functions to load user config files. These are called by the main MultiQC script.
# Note that config files are loaded in a specific order and values can overwrite each other.
def mqc_load_userconfig(paths=()):
    """Overwrite config defaults with user config files"""

    # Load and parse installation config file if we find it
    mqc_load_config(os.path.join(os.path.dirname(MULTIQC_DIR), "multiqc_config.yaml"))

    # Load and parse a config file in $XDG_CONFIG_HOME
    # Ref: https://specifications.freedesktop.org/basedir-spec/basedir-spec-latest.html
    mqc_load_config(
        os.path.join(os.environ.get("XDG_CONFIG_HOME", os.path.expanduser("~/.config")), "multiqc_config.yaml")
    )

    # Load and parse a user config file if we find it
    mqc_load_config(os.path.expanduser("~/.multiqc_config.yaml"))

    # Load and parse a config file path set in an ENV variable if we find it
    if os.environ.get("MULTIQC_CONFIG_PATH") is not None:
        mqc_load_config(os.environ.get("MULTIQC_CONFIG_PATH"))

    # Load separate config entries from MULTIQC_* environment variables
    mqc_add_config(mqc_env_vars_config())

    # Load and parse a config file in this working directory if we find it
    mqc_load_config("multiqc_config.yaml")

    # Custom command line config
    for p in paths:
        mqc_load_config(p)


def mqc_load_config(yaml_config_path: str):
    """Load and parse a config file if we find it"""
    if not os.path.isfile(yaml_config_path) and os.path.isfile(yaml_config_path.replace(".yaml", ".yml")):
        yaml_config_path = yaml_config_path.replace(".yaml", ".yml")

    if os.path.isfile(yaml_config_path):
        try:
            # pyaml_env allows referencing environment variables in YAML for default values
            new_config = pyaml_env.parse_config(yaml_config_path)
            logger.debug(f"Loading config settings from: {yaml_config_path}")
            mqc_add_config(new_config, yaml_config_path)
        except (IOError, AttributeError) as e:
            logger.debug(f"Config error: {e}")
        except yaml.scanner.ScannerError as e:
            logger.error(f"Error parsing config YAML: {e}")
            sys.exit(1)


def mqc_cl_config(cl_config):
    for clc_str in cl_config:
        try:
            parsed_clc = yaml.safe_load(clc_str)
            # something:var fails as it needs a space. Fix this (a common mistake)
            if isinstance(parsed_clc, str) and ":" in clc_str:
                clc_str = ": ".join(clc_str.split(":"))
                parsed_clc = yaml.safe_load(clc_str)
            assert isinstance(parsed_clc, dict)
        except yaml.scanner.ScannerError as e:
            logger.error(f"Could not parse command line config: {clc_str}\n{e}")
        except AssertionError:
            logger.error(f"Could not parse command line config: {clc_str}")
        else:
            logger.debug(f"Found command line config: {parsed_clc}")
            mqc_add_config(parsed_clc)


def mqc_env_vars_config() -> Dict:
    """
    Check MULTIQC_* environment variables and set to corresponding config values if they are of scalar types.
    """
    RESERVED_NAMES = {"MULTIQC_CONFIG_PATH"}
    PREFIX = "MULTIQC_"  # Prefix for environment variables
    env_config = {}
    for k, v in os.environ.items():
        if k.startswith(PREFIX) and k not in RESERVED_NAMES:
            conf_key = k[len(PREFIX) :].lower()
            if conf_key not in globals():
                continue
            if isinstance(globals()[conf_key], bool):
                try:
                    v = strtobool(v)
                except ValueError:
                    logger.warning(f"Could not parse a boolean value from the environment variable ${k}={v}")
                    continue
            elif isinstance(globals()[conf_key], int):
                try:
                    v = int(v)
                except ValueError:
                    logger.warning(f"Could not parse a int value from the environment variable ${k}={v}")
                    continue
            elif isinstance(globals()[conf_key], float):
                try:
                    v = float(v)
                except ValueError:
                    logger.warning(f"Could not parse a float value from the environment variable ${k}={v}")
                    continue
            elif not isinstance(globals()[conf_key], str) and globals()[conf_key] is not None:
                logger.warning(
                    f"Can only set scalar config entries (str, int, float, bool) with environment variable, "
                    f"but config.{conf_key} expects a type '{type(globals()[conf_key]).__name__}'. Ignoring ${k}"
                )
                continue
            env_config[conf_key] = v
            logger.debug(f"Setting config.{conf_key} from the environment variable ${k}")
    return env_config


def mqc_add_config(conf: Dict, conf_path=None):
    """Add to the global config with given MultiQC config dict"""
    global custom_css_files, fn_clean_exts, fn_clean_trim
    log_new_config = {}
    log_filename_patterns = []
    log_filename_clean_extensions = []
    log_filename_clean_trimmings = []
    for c, v in conf.items():
        if c == "sp":
            # Merge filename patterns instead of replacing
            sp.update(v)
            log_filename_patterns.append(v)
        elif c == "extra_fn_clean_exts":
            # Prepend to filename cleaning patterns instead of replacing
            fn_clean_exts[0:0] = v
            log_filename_clean_extensions.append(v)
        elif c == "extra_fn_clean_trim":
            # Prepend to filename cleaning patterns instead of replacing
            fn_clean_trim[0:0] = v
            log_filename_clean_trimmings.append(v)
        elif c in ["custom_logo"] and v:
            # Resolve file paths - absolute or cwd, or relative to config file
            fpath = v
            if os.path.exists(v):
                fpath = os.path.abspath(v)
            elif conf_path is not None and os.path.exists(os.path.join(os.path.dirname(conf_path), v)):
                fpath = os.path.abspath(os.path.join(os.path.dirname(conf_path), v))
            else:
                logger.error(f"Config '{c}' path not found, skipping ({fpath})")
                continue
            log_new_config[c] = fpath
            update({c: fpath})
        elif c == "custom_css_files":
            for fpath in v:
                if os.path.exists(fpath):
                    fpath = os.path.abspath(fpath)
                elif conf_path is not None and os.path.exists(os.path.join(os.path.dirname(conf_path), fpath)):
                    fpath = os.path.abspath(os.path.join(os.path.dirname(conf_path), fpath))
                else:
                    logger.error(f"CSS path '{c}' path not found, skipping ({fpath})")
                    continue
                logger.debug(f"Adding css file '{c}': {fpath}")
                if not custom_css_files:
                    custom_css_files = []
                custom_css_files.append(fpath)
        else:
            log_new_config[c] = v
            update({c: v})
    if len(log_new_config) > 0:
        logger.debug(f"New config: {log_new_config}")
    if len(log_filename_patterns) > 0:
        logger.debug(f"Added to filename patterns: {log_filename_patterns}")
    if len(log_filename_clean_extensions) > 0:
        logger.debug(f"Added to filename clean extensions: {log_filename_clean_extensions}")
    if len(log_filename_clean_trimmings) > 0:
        logger.debug(f"Added to filename clean trimmings: {log_filename_clean_trimmings}")


#### Function to load file containing a list of alternative sample-name swaps
# Essentially a fancy way of loading stuff into the sample_names_rename config var
# As such, can also be done directly using a config file
def load_sample_names(snames_file):
    global sample_names_rename_buttons, sample_names_rename
    num_cols = None
    try:
        with open(snames_file) as f:
            logger.debug(f"Loading sample renaming config settings from: {snames_file}")
            for line in f:
                s = line.strip().split("\t")
                if len(s) > 1:
                    # Check that we have consistent numbers of columns
                    if num_cols is None:
                        num_cols = len(s)
                    elif num_cols != len(s):
                        logger.warning(
                            "Inconsistent number of columns found in sample names file (skipping line): '{}'".format(
                                line.strip()
                            )
                        )
                    # Parse the line
                    if len(sample_names_rename_buttons) == 0:
                        sample_names_rename_buttons = s
                    else:
                        sample_names_rename.append(s)
                elif len(line.strip()) > 0:
                    logger.warning(f"Sample names file line did not have columns (must use tabs): {line.strip()}")
    except (IOError, AttributeError) as e:
        logger.error(f"Error loading sample names file: {e}")
    logger.debug(f"Found {len(sample_names_rename_buttons)} sample renaming patterns")


def load_replace_names(rnames_file):
    global sample_names_replace
    try:
        with open(rnames_file) as f:
            logger.debug(f"Loading sample replace config settings from: {rnames_file}")
            for line in f:
                s = line.strip().split("\t")
                if len(s) == 2:
                    sample_names_replace[s[0]] = s[1]
    except (IOError, AttributeError) as e:
        logger.error(f"Error loading sample names replacement file: {e}")
    logger.debug(f"Found {len(sample_names_replace)} sample replacing patterns")


def load_show_hide(sh_file):
    global show_hide_buttons, show_hide_patterns, show_hide_mode, show_hide_regex
    if sh_file:
        try:
            with open(sh_file, "r") as f:
                logger.debug(f"Loading sample renaming config settings from: {sh_file}")
                for line in f:
                    s = line.strip().split("\t")
                    if len(s) >= 3 and s[1] in ["show", "hide", "show_re", "hide_re"]:
                        show_hide_buttons.append(s[0])
                        show_hide_mode.append(s[1])
                        show_hide_patterns.append(s[2:])
                        show_hide_regex.append(s[1] not in ["show", "hide"])  # flag whether or not regex is turned on
        except AttributeError as e:
            logger.error(f"Error loading show patterns file: {e}")

    # Prepend a "Show all" button if we have anything
    # Do this outside of the file load block in case it was set in the config
    if len(show_hide_buttons) > 0:
        logger.debug(f"Found {len(show_hide_buttons)} show/hide patterns")
        show_hide_buttons.insert(0, "Show all")
        show_hide_mode.insert(0, "hide")
        show_hide_patterns.insert(0, [])
        show_hide_regex.insert(0, False)


# Keep track of all changes to the config
nondefault_config = dict()


def update(u):
    update_dict(nondefault_config, u)
    return update_dict(globals(), u)


def update_dict(target, source, none_only=False):
    """Recursively updates nested dict d from nested dict u"""
    for key, val in source.items():
        if isinstance(val, dict):
            target[key] = update_dict(target.get(key, {}), val)
        else:
            if not none_only or target.get(key) is None:
                if isinstance(val, list):
                    target[key] = val.copy()
                else:
                    target[key] = val
    return target
