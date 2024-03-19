""" MultiQC Utility functions, used in a variety of places. """

import io
import json
import logging
from collections import defaultdict, OrderedDict

import math
import os
import shutil
import sys
import time
import datetime
from typing import Dict, List, Union

import yaml

from . import config

log = logging.getLogger(__name__)


def robust_rmtree(path, logger=None, max_retries=10):
    """Robustly tries to delete paths.
    Retries several times (with increasing delays) if an OSError
    occurs.  If the final attempt fails, the Exception is propagated
    to the caller.
    """

    for i in range(max_retries):
        try:
            shutil.rmtree(path)
            return
        except OSError:
            if logger:
                logger.info(f"Unable to remove path: {path}")
                logger.info(f"Retrying after {i**2} seconds")
            else:
                print(f"Unable to remove path: {path}", file=sys.stderr)
                print(f"Retrying after {i**2} seconds", file=sys.stderr)
            time.sleep(i**2)

    # Final attempt, pass any Exceptions up to caller.
    shutil.rmtree(path)


def write_data_file(
    data: Union[Dict[str, Union[Dict, List]], List[Dict]],
    fn: str,
    sort_cols=False,
    data_format=None,
):
    """
    Write a data file to the report directory. Will not do anything
    if config.data_dir is not set.
    :param: data - either: a 2D dict, first key sample name (row header),
        second key field (column header); a list of dicts; or a list of lists
    :param: fn - Desired filename. Directory will be prepended automatically.
    :param: sort_cols - Sort columns alphabetically
    :param: data_format - Output format. Defaults to config.data_format (usually tsv)
    :return: None
    """

    if config.data_dir is None:
        return

    # Get data format from config
    if data_format is None:
        data_format = config.data_format

    body = None
    # Some metrics can't be coerced to tab-separated output, test and handle exceptions
    if data_format in ["tsv", "csv"]:
        sep = "\t" if data_format == "tsv" else ","
        # Attempt to reshape data to tsv
        # noinspection PyBroadException
        try:
            # Get all headers from the data, except if data is a dictionary (i.e. has >1 dimensions)
            headers = []
            rows = []

            for d in data.values() if isinstance(data, dict) else data:
                if not d or (isinstance(d, list) and isinstance(d[0], dict)):
                    continue
                if isinstance(d, dict):
                    for h in d.keys():
                        if h not in headers:
                            headers.append(h)
            if headers:
                if sort_cols:
                    headers = sorted(headers)
                headers_str = [str(item) for item in headers]
                if isinstance(data, dict):
                    # Add Sample header as a first element
                    headers_str.insert(0, "Sample")
                rows.append(sep.join(headers_str))

            # The rest of the rows
            for key, d in sorted(data.items()) if isinstance(data, dict) else enumerate(data):
                # Make a list starting with the sample name, then each field in order of the header cols
                if headers:
                    line = [str(d.get(h, "")) for h in headers]
                else:
                    line = [
                        str(item)
                        for item in (d.values() if isinstance(d, dict) else (d if isinstance(d, list) else [d]))
                    ]
                if isinstance(data, dict):
                    # Add Sample header as a first element
                    line.insert(0, str(key))
                rows.append(sep.join(line))
            body = "\n".join(rows)

        except Exception as e:
            if config.development:
                raise
            data_format = "yaml"
            log.debug(f"{fn} could not be saved as tsv/csv, falling back to YAML. {e}")

    # Add relevant file extension to filename, save file.
    fn = f"{fn}.{config.data_format_extensions[data_format]}"
    fpath = os.path.join(config.data_dir, fn)
    with io.open(fpath, "w", encoding="utf-8") as f:
        if data_format == "json":
            jsonstr = dump_json(data, indent=4, ensure_ascii=False)
            print(jsonstr.encode("utf-8", "ignore").decode("utf-8"), file=f)
        elif data_format == "yaml":
            yaml.dump(data, f, default_flow_style=False)
        elif body:
            # Default - tab separated output
            print(body.encode("utf-8", "ignore").decode("utf-8"), file=f)
    log.debug(f"Wrote data file {fn}")


def force_term_colors():
    """
    Check if any environment variables are set to force Rich to use coloured output
    """
    if os.getenv("GITHUB_ACTIONS") or os.getenv("FORCE_COLOR") or os.getenv("PY_COLORS"):
        return True
    return None


def strtobool(val) -> bool:
    """
    Replaces deprecated https://docs.python.org/3.9/distutils/apiref.html#distutils.util.strtobool
    The deprecation recommendation is to re-implement the function https://peps.python.org/pep-0632/

    ------------------------------------------------------------

    Convert a string representation of truth to true (1) or false (0).

    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.
    """
    val_str = str(val).lower()
    if val_str in ("y", "yes", "t", "true", "on", "1"):
        return True
    elif val_str in ("n", "no", "f", "false", "off", "0"):
        return False
    else:
        raise ValueError(f"invalid truth value {val!r}")


def choose_emoji():
    """Choose an emoji to use in the report header."""
    # NB: We haven't parsed the config yet, so can't disable via config
    today = datetime.date.today()
    emojis = {
        "bottle_with_popping_cork": (1, 1, 1, 5),  # New Year's Day
        "rose": (2, 14, 0, 0),  # Valentine's Day
        "four_leaf_clover": (3, 17, 0, 0),  # St Patrick's Day
        "globe_showing_asia-australia": (4, 22, 0, 0),  # Earth Day
        "jack-o-lantern": (10, 31, 5, 0),  # Halloween
        "santa": (12, 25, 0, 0),  # Christmas Day
        "christmas_tree": (12, 25, 7, 7),  # Christmas
    }
    for emoji, (month, day, days_before, days_after) in emojis.items():
        special_date = datetime.date(today.year, month, day)
        date_range_start = special_date - datetime.timedelta(days=days_before)
        date_range_end = special_date + datetime.timedelta(days=days_after)
        if date_range_start <= today <= date_range_end:
            return emoji
    return "mag"


def dump_json(data, **kwargs):
    """
    Recursively replace non-JSON-conforming NaNs and lambdas with None.
    Note that a custom JSONEncoder would have worked for lambdas, but not for NaNs: https://stackoverflow.com/a/28640141
    """

    # Recursively replace NaNs with None
    def replace_nan(obj):
        if isinstance(obj, dict):
            return {k: replace_nan(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [replace_nan(v) for v in obj]
        elif isinstance(obj, set):
            return {replace_nan(v) for v in obj}
        elif callable(obj):
            return None
        elif isinstance(obj, float) and math.isnan(obj):
            return None
        return obj

    return json.dumps(replace_nan(data), **kwargs)


def multiqc_dump_json(report):
    """
    Export the parsed data in memory to a JSON file.
    Used for MegaQC and other data export.
    WARNING: May be depreciated and removed in future versions.
    """
    exported_data = dict()
    export_vars = {
        "report": [
            "data_sources",
            "general_stats_data",
            "general_stats_headers",
            "multiqc_command",
            "plot_data",
            "saved_raw_data",
        ],
        "config": [
            "analysis_dir",
            "creation_date",
            "git_hash",
            "intro_text",
            "report_comment",
            "report_header_info",
            "script_path",
            "short_version",
            "subtitle",
            "title",
            "version",
            "output_dir",
        ],
    }
    for s in export_vars:
        for k in export_vars[s]:
            try:
                d = None
                if s == "config":
                    d = {f"{s}_{k}": getattr(config, k)}
                elif s == "report":
                    d = {f"{s}_{k}": getattr(report, k)}
                if d:
                    dump_json(d, ensure_ascii=False)  # Test that exporting to JSON works
                    exported_data.update(d)
            except (TypeError, KeyError, AttributeError) as e:
                log.warning(f"Couldn't export data key '{s}.{k}': {e}")
        # Get the absolute paths of analysis directories
        exported_data["config_analysis_dir_abs"] = list()
        for d in exported_data.get("config_analysis_dir", []):
            try:
                exported_data["config_analysis_dir_abs"].append(os.path.abspath(d))
            except Exception:
                pass
    return exported_data


def replace_defaultdicts(data):
    """
    Recursively replace dict-likes as dicts for nice yaml representation.
    """

    def _replace(obj):
        if isinstance(obj, (defaultdict, OrderedDict, dict)):
            return {k: _replace(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [_replace(v) for v in obj]
        elif isinstance(obj, set):
            return {_replace(v) for v in obj}
        return obj

    return _replace(data)
