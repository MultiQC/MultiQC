#!/usr/bin/env python

""" MultiQC Utility functions, used in a variety of places. """

from __future__ import print_function
from collections import OrderedDict
import io
import json
from logging import NullHandler
import os
import yaml
import time
import shutil
import sys

from multiqc import config


def tab_separated_output(data, sort_cols):
    """Subfunction of write_data_file() below
    Some metrics, in particular those of the custom content module
    can't be written cleanly to tab-separated output. Now failure
    can trigger the use of a more suitable export format.
    """
    try:
        # Convert keys to strings
        data = {str(k): v for k, v in data.items()}
        # Get all headers
        h = ["Sample"]
        for sn in sorted(data.keys()):
            for k in data[sn].keys():
                if type(data[sn][k]) is not dict and k not in h:
                    h.append(str(k))
        if sort_cols:
            h = sorted(h)

        # Get the rows
        rows = ["\t".join(h)]
        for sn in sorted(data.keys()):
            # Make a list starting with the sample name, then each field in order of the header cols
            l = [str(sn)] + [str(data[sn].get(k, "")) for k in h[1:]]
            rows.append("\t".join(l))

        body = "\n".join(rows)

        return body, True
    except:
        return None, False


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
                logger.info("Unable to remove path: {}".format(path))
                logger.info("Retrying after {} seconds".format(i ** 2))
            else:
                print("Unable to remove path: {}".format(path), file=sys.stderr)
                print("Retrying after {} seconds".format(i ** 2), file=sys.stderr)
            time.sleep(i ** 2)

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

        # Initialise variables
        fallback = False
        success = False

        # Get data format from config
        if data_format is None:
            data_format = config.data_format

        # JSON encoder class to handle lambda functions
        class MQCJSONEncoder(json.JSONEncoder):
            def default(self, obj):
                if callable(obj):
                    try:
                        return obj(1)
                    except:
                        return None
                return json.JSONEncoder.default(self, obj)

        # Some metrics can't be coerced to tab-separated output, test and handle exceptions
        if data_format not in ["json", "yaml"]:

            # attempt to reshape data to tsv
            data_tsv, success = tab_separated_output(data, sort_cols=sort_cols)

            if not success:
                # Unsuccessful TSV reshape: Try to remove malformed data and export remainders.
                aberrant_data = list()
                regular_data = list()

                try:
                    for _, key in zip(range(len(data)), data):
                        _, success = tab_separated_output([data[key]], sort_cols=sort_cols)
                        if success:
                            regular_data.append(data[key])
                        else:
                            aberrant_data.append(data[key])

                    if regular_data:
                        # Final, now hopefully working reshape
                        data_tsv, success = tab_separated_output(regular_data, sort_cols=sort_cols)

                        if success:
                            # Separate export for what didn't work.
                            for key in aberrant_data:
                                adf = str(fn + "_" + str(key) + config.data_format_extensions["yaml"])
                                with io.open(os.path.join(config.data_dir, adf), "w", encoding="utf-8") as g:
                                    yaml.dump(aberrant_data[key], g, default_flow_style=False)
                                    config.logger.info(
                                        "Some metrics could not be coerced into tab-separated output and were instead written as YAML to "
                                        + adf
                                    )
                        else:  # Reshaping remainder failed.
                            fallback = True
                    else:  # No regular data: Rather export everything as YAML.
                        fallback = True
                except:
                    # Subsetting doesn't work, e.g. because key is a list for scatterplots
                    fallback = True

        if fallback:
            data_format = "yaml"
            config.logger.warning(
                "Metrics of "
                + fn
                + " could not be saved as tab-separated output. Choose JSON or YAML output instead. Falling back to YAML."
            )

        # Add relevant file extension to filename, save file.
        fn = "{}.{}".format(fn, config.data_format_extensions[data_format])
        with io.open(os.path.join(config.data_dir, fn), "w", encoding="utf-8") as f:
            if data_format == "json":
                jsonstr = json.dumps(data, indent=4, cls=MQCJSONEncoder, ensure_ascii=False)
                print(jsonstr.encode("utf-8", "ignore").decode("utf-8"), file=f)
            elif data_format == "yaml":
                yaml.dump(data, f, default_flow_style=False)
            else:
                # Default - tab separated output
                print(data_tsv.encode("utf-8", "ignore").decode("utf-8"), file=f)


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
        print(" - {}:".format(t))
        for ttgs in avail_tags[t]:
            print("   - {}".format(ttgs))
    ctx.exit()
