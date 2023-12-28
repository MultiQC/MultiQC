#!/usr/bin/env python

""" MultiQC Utility functions, used in a variety of places. """


import io
import json
import os
import shutil
import sys
import time

import yaml

from . import config


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


def write_data_file(data, fn, sort_cols=False, data_format=None):
    """Write a data file to the report directory. Will not do anything
    if config.data_dir is not set.
    :param: data - a 2D dict, first key sample name (row header),
            second key field (column header).
    :param: fn - Desired filename. Directory will be prepended automatically.
    :param: sort_cols - Sort columns alphabetically
    :param: data_format - Output format. Defaults to config.data_format (usually tsv)
    :return: None"""

    if config.data_dir is not None:
        # Get data format from config
        if data_format is None:
            data_format = config.data_format

        # JSON encoder class to handle lambda functions
        class MQCJSONEncoder(json.JSONEncoder):
            def default(self, obj):
                if callable(obj):
                    try:
                        return obj(1)
                    except Exception:
                        return None
                return json.JSONEncoder.default(self, obj)

        # Some metrics can't be coerced to tab-separated output, test and handle exceptions
        if data_format not in ["json", "yaml"]:
            # attempt to reshape data to tsv
            try:
                # Get all headers from the data, except if data is a dictionary (i.e. has >1 dimensions)
                headers = []
                for d in data.values():
                    if not d or (isinstance(d, list) and isinstance(d[0], dict)):
                        continue
                    for h in d.keys():
                        if h not in headers:
                            headers.append(h)
                if sort_cols:
                    headers = sorted(headers)

                headers_str = [str(item) for item in headers]
                # Add Sample header in to first element
                headers_str.insert(0, "Sample")

                # Get the rows
                rows = ["\t".join(headers_str)]
                for sn in sorted(data.keys()):
                    # Make a list starting with the sample name, then each field in order of the header cols
                    line = [str(sn)] + [str(data[sn].get(h, "")) for h in headers]
                    rows.append("\t".join(line))

                body = "\n".join(rows)

            except Exception:
                data_format = "yaml"
                config.logger.debug(f"{fn} could not be saved as tsv/csv. Falling back to YAML.")

        # Add relevant file extension to filename, save file.
        fn = f"{fn}.{config.data_format_extensions[data_format]}"
        with io.open(os.path.join(config.data_dir, fn), "w", encoding="utf-8") as f:
            if data_format == "json":
                jsonstr = json.dumps(data, indent=4, cls=MQCJSONEncoder, ensure_ascii=False)
                print(jsonstr.encode("utf-8", "ignore").decode("utf-8"), file=f)
            elif data_format == "yaml":
                yaml.dump(data, f, default_flow_style=False)
            else:
                # Default - tab separated output
                print(body.encode("utf-8", "ignore").decode("utf-8"), file=f)


def view_all_tags(ctx, param, value):
    """List available tags and associated modules
    Called by eager click option: --view-tags
    """
    # To make sure this function executed only when the flag was called
    if not value or ctx.resilient_parsing:
        return
    avail_tags = dict()
    print("\nMultiQC Available module tag groups:\n")
    for mod_dict in filter(lambda mod: isinstance(mod, dict), config.module_order):
        mod_key, mod_val = list(mod_dict.items())[0]
        tags = list(mod_val.get("module_tag", []))
        for t in tags:
            if t not in avail_tags:
                avail_tags[t] = []
            avail_tags[t].append(mod_key)
    for t in sorted(avail_tags.keys(), key=lambda s: s.lower()):
        print(f" - {t}:")
        for ttgs in avail_tags[t]:
            print(f"   - {ttgs}")
    ctx.exit()


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
    val = val.lower()
    if val in ("y", "yes", "t", "true", "on", "1"):
        return True
    elif val in ("n", "no", "f", "false", "off", "0"):
        return False
    else:
        raise ValueError(f"invalid truth value {val!r}")
