""" MultiQC datatable class, used by tables and violin plots """

import math

import logging
import random
import re
import string
from collections import defaultdict
from typing import List, Tuple, Dict, Optional, Union

from multiqc.utils import config, report

logger = logging.getLogger(__name__)


class DataTable:
    """
    Data table class. Prepares and holds data and configuration
    for either a table or a violin plot.
    """

    def __init__(
        self,
        data: Union[List[Dict], Dict],
        headers: Optional[Union[List, Dict]] = None,
        pconfig: Optional[Dict] = None,
    ):
        """Prepare data for use in a table or plot"""
        self.headers_in_order = defaultdict(list)
        self.data: Dict = {}
        self.headers: Optional[List] = None
        self.pconfig: Optional[Dict] = None

        if headers is None:
            headers = []
        if pconfig is None:
            pconfig = {}

        # Given one dataset - turn it into a list
        if not isinstance(data, list):
            data = [data]
        if not isinstance(headers, list):
            headers = [headers]

        if pconfig and "id" in pconfig:
            self.id = pconfig.pop("id")
        else:
            if config.strict:
                errmsg = f"LINT: 'id' is missing from plot pconfig: {pconfig}"
                logger.error(errmsg)
                report.lint_errors.append(errmsg)
            self.id = report.save_htmlid(f"table-{''.join(random.sample(string.ascii_lowercase, 4))}")

        self._build(data, headers, pconfig)

    def _build(
        self,
        data: List[Dict],
        headers: List[Dict],
        pconfig: Dict,
    ):
        # Allow user to overwrite any given config for this plot
        if self.id in config.custom_plot_config:
            for k, v in config.custom_plot_config[self.id].items():
                pconfig[k] = v

        SECTION_COLORS = [
            "55,126,184",  # Blue
            "77,175,74",  # Green
            "152,78,163",  # Purple
            "255,127,0",  # Orange
            "228,26,28",  # Red
            "179,179,50",  # Olive
            "166,86,40",  # Brown
            "247,129,191",  # Pink
            "153,153,153",  # Grey
        ]

        # Go through each table section
        for idx, d in enumerate(data):
            # Get the header keys
            try:
                keys = headers[idx].keys()
                assert len(keys) > 0
            except (IndexError, AttributeError, AssertionError):
                pconfig["only_defined_headers"] = False
                keys = list()

            # Add header keys from the data
            if pconfig.get("only_defined_headers", True) is False:
                # Get the keys from the data
                keys = list()
                for samp in d.values():
                    for k in samp.keys():
                        if k not in keys:
                            keys.append(k)

                # If we don't have a headers dict for this data set yet, create one
                try:
                    headers[idx]
                except IndexError:
                    headers.append(dict())
                else:
                    # Convert the existing headers into a dict (e.g. if parsed from a config)
                    od_tuples = [(key, headers[idx][key]) for key in headers[idx].keys()]
                    headers[idx] = dict(od_tuples)

                # Create empty header configs for each new data key
                for k in keys:
                    if k not in headers[idx]:
                        headers[idx][k] = {}

            # Ensure that keys are strings, not numeric
            keys = [str(k) for k in keys]
            for k in list(headers[idx].keys()):
                headers[idx][str(k)] = headers[idx].pop(k)
            # Ensure that all sample names are strings as well
            cdata = dict()
            for k, v in data[idx].items():
                cdata[str(k)] = v
            data[idx] = cdata
            for s_name in data[idx].keys():
                for k in list(data[idx][s_name].keys()):
                    data[idx][s_name][str(k)] = data[idx][s_name].pop(k)

            # Check that we have some data in each column
            empties = list()
            for k in keys:
                n = 0
                for samp in d.values():
                    if k in samp:
                        n += 1
                if n == 0:
                    empties.append(k)
            for k in empties:
                keys = [j for j in keys if j != k]
                del headers[idx][k]

            for k in keys:
                # Unique id to avoid overwriting by other datasets
                if "rid" not in headers[idx][k]:
                    headers[idx][k]["rid"] = report.save_htmlid(re.sub(r"\W+", "_", k).strip().strip("_"))

                # Applying defaults presets for data keys if shared_key is set to base_count or read_count
                shared_key = headers[idx][k].get("shared_key", None)
                shared_key_suffix = None
                if shared_key in ["read_count", "long_read_count", "base_count"]:
                    if shared_key == "read_count" and config.read_count_prefix:
                        multiplier = config.read_count_multiplier
                        shared_key_suffix = config.read_count_prefix
                    elif shared_key == "long_read_count" and config.long_read_count_prefix:
                        multiplier = config.long_read_count_multiplier
                        shared_key_suffix = config.long_read_count_prefix
                    elif shared_key == "base_count" and config.base_count_prefix:
                        multiplier = config.base_count_multiplier
                        shared_key_suffix = config.base_count_prefix
                    else:
                        multiplier = 1
                    if "modify" not in headers[idx][k]:
                        headers[idx][k]["modify"] = lambda x: x * multiplier
                    if "min" not in headers[idx][k] is None:
                        headers[idx][k]["min"] = 0
                    if "format" not in headers[idx][k] is None:
                        if multiplier == 1:
                            headers[idx][k]["format"] = "{:,d}"
                if "suffix" not in headers[idx][k] and shared_key_suffix is not None:
                    headers[idx][k]["suffix"] = " " + shared_key_suffix

                # Use defaults / data keys if headers not given
                headers[idx][k]["namespace"] = headers[idx][k].get("namespace", pconfig.get("namespace", ""))
                headers[idx][k]["title"] = headers[idx][k].get("title", k)
                headers[idx][k]["description"] = headers[idx][k].get("description", headers[idx][k]["title"])
                headers[idx][k]["scale"] = headers[idx][k].get("scale", pconfig.get("scale", "GnBu"))
                headers[idx][k]["format"] = headers[idx][k].get("format", pconfig.get("format", "{:,.1f}"))
                headers[idx][k]["colour"] = headers[idx][k].get("colour", pconfig.get("colour", None))
                headers[idx][k]["hidden"] = headers[idx][k].get("hidden", pconfig.get("hidden", None))
                headers[idx][k]["max"] = headers[idx][k].get("max", pconfig.get("max", None))
                headers[idx][k]["min"] = headers[idx][k].get("min", pconfig.get("min", None))
                headers[idx][k]["ceiling"] = headers[idx][k].get("ceiling", pconfig.get("ceiling", None))
                headers[idx][k]["floor"] = headers[idx][k].get("floor", pconfig.get("floor", None))
                headers[idx][k]["minrange"] = headers[idx][k].get(
                    "minrange", pconfig.get("minrange", headers[idx][k].get("minRange", pconfig.get("minRange", None)))
                )
                headers[idx][k]["shared_key"] = headers[idx][k].get("shared_key", pconfig.get("shared_key", None))
                headers[idx][k]["modify"] = headers[idx][k].get("modify", pconfig.get("modify", None))
                headers[idx][k]["placement"] = float(headers[idx][k].get("placement", 1000))

                if headers[idx][k]["colour"] is None:
                    cidx = idx
                    while cidx >= len(SECTION_COLORS):
                        cidx -= len(SECTION_COLORS)
                    headers[idx][k]["colour"] = SECTION_COLORS[cidx]

                # Overwrite (2nd time) any given config with table-level user config
                # This is to override column-specific values set by modules
                if self.id in config.custom_plot_config:
                    for cpc_k, cpc_v in config.custom_plot_config[self.id].items():
                        headers[idx][k][cpc_k] = cpc_v

                # Overwrite "name" if set in user config
                # Key can be a column ID, a table ID, or a namespace in the general stats table.
                for key, val in config.table_columns_name.items():
                    key = key.lower()
                    # Case-insensitive check if the outer key is a table ID or a namespace.
                    if key in [self.id.lower(), headers[idx][k]["namespace"].lower()] and isinstance(val, dict):
                        # Assume a dict of specific column IDs
                        for key2, new_title in val.items():
                            key2 = key2.lower()
                            if key2 in [k.lower(), headers[idx][k]["title"].lower()]:
                                headers[idx][k]["title"] = new_title

                    # Case-insensitive check if the outer key is a column ID
                    elif key in [k.lower(), headers[idx][k]["title"].lower()] and isinstance(val, str):
                        headers[idx][k]["title"] = val

                # Overwrite "hidden" if set in user config
                # Key can be a column ID, a table ID, or a namespace in the general stats table.
                for key, val in config.table_columns_visible.items():
                    key = key.lower()
                    # Case-insensitive check if the outer key is a table ID or a namespace.
                    if key in [self.id.lower(), headers[idx][k]["namespace"].lower()]:
                        # First - if config value is a bool, set all module columns to that value
                        if isinstance(val, bool):
                            # Config has True = visible, False = Hidden. Here we're setting "hidden" which is inverse
                            headers[idx][k]["hidden"] = not val

                        # Not a bool, assume a dict of specific column IDs
                        elif isinstance(val, dict):
                            for key2, visible in val.items():
                                key2 = key2.lower()
                                if key2 in [k.lower(), headers[idx][k]["title"].lower()] and isinstance(visible, bool):
                                    # Config has True = visible, False = Hidden. Here we're setting "hidden" which is inverse
                                    headers[idx][k]["hidden"] = not visible

                    # Case-insensitive check if the outer key is a column ID
                    elif key in [k.lower(), headers[idx][k]["title"].lower()] and isinstance(val, bool):
                        # Config has True = visible, False = Hidden. Here we're setting "hidden" which is inverse
                        headers[idx][k]["hidden"] = not val

                # Also overwrite placement if set in config
                try:
                    headers[idx][k]["placement"] = float(
                        config.table_columns_placement[headers[idx][k]["namespace"]][k]
                    )
                except (KeyError, ValueError):
                    try:
                        headers[idx][k]["placement"] = float(config.table_columns_placement[self.id][k])
                    except (KeyError, ValueError):
                        pass

                # Overwrite any header config if set in config
                for custom_k, custom_v in config.custom_table_header_config.get(self.id, {}).get(k, {}).items():
                    headers[idx][k][custom_k] = custom_v

                # Work out max and min value if not given
                setdmax = False
                setdmin = False
                try:
                    headers[idx][k]["dmax"] = float(headers[idx][k]["max"])
                except Exception:
                    headers[idx][k]["dmax"] = 0
                    setdmax = True

                try:
                    headers[idx][k]["dmin"] = float(headers[idx][k]["min"])
                except Exception:
                    headers[idx][k]["dmin"] = 0
                    setdmin = True

                # Figure out the min / max if not supplied
                if setdmax or setdmin:
                    for s_name, samp in data[idx].items():
                        try:
                            val = float(samp[k])
                            if callable(headers[idx][k]["modify"]):
                                val = float(headers[idx][k]["modify"](val))
                            if math.isfinite(val) and not math.isnan(val):
                                if setdmax:
                                    headers[idx][k]["dmax"] = max(headers[idx][k]["dmax"], val)
                                if setdmin:
                                    headers[idx][k]["dmin"] = min(headers[idx][k]["dmin"], val)
                        except (ValueError, TypeError):
                            val = samp[k]  # couldn't convert to float - keep as a string
                        except KeyError:
                            pass  # missing data - skip
                    # Limit auto-generated scales with floor, ceiling and minrange.
                    if headers[idx][k]["ceiling"] is not None and headers[idx][k]["max"] is None:
                        headers[idx][k]["dmax"] = min(headers[idx][k]["dmax"], float(headers[idx][k]["ceiling"]))
                    if headers[idx][k]["floor"] is not None and headers[idx][k]["min"] is None:
                        headers[idx][k]["dmin"] = max(headers[idx][k]["dmin"], float(headers[idx][k]["floor"]))
                    if headers[idx][k]["minrange"] is not None:
                        drange = headers[idx][k]["dmax"] - headers[idx][k]["dmin"]
                        if drange < float(headers[idx][k]["minrange"]):
                            headers[idx][k]["dmax"] = headers[idx][k]["dmin"] + float(headers[idx][k]["minrange"])

        # Collect settings for shared keys
        shared_keys = defaultdict(lambda: dict())
        for idx, hs in enumerate(headers):
            for k in hs.keys():
                sk = headers[idx][k]["shared_key"]
                if sk is not None:
                    shared_keys[sk]["dmax"] = max(
                        headers[idx][k]["dmax"], shared_keys[sk].get("dmax", headers[idx][k]["dmax"])
                    )
                    shared_keys[sk]["dmin"] = min(
                        headers[idx][k]["dmin"], shared_keys[sk].get("dmin", headers[idx][k]["dmin"])
                    )

        # Overwrite shared key settings and at the same time assign to buckets for sorting
        # So the final ordering is:
        #   placement > section > explicit_ordering
        # Of course, the user can shuffle these manually.
        for idx, hs in enumerate(headers):
            for k in hs.keys():
                sk = headers[idx][k]["shared_key"]
                if sk is not None:
                    headers[idx][k]["dmax"] = shared_keys[sk]["dmax"]
                    headers[idx][k]["dmin"] = shared_keys[sk]["dmin"]

                self.headers_in_order[headers[idx][k]["placement"]].append((idx, k))

        # Skip any data that is not used in the table
        # Would be ignored for making the table anyway, but can affect whether a violin plot is used
        for idx, d in enumerate(data):
            for s_name in list(d.keys()):
                if not any(h in data[idx][s_name].keys() for h in headers[idx]):
                    del data[idx][s_name]

        # Assign to class
        self.data = data
        self.headers = headers
        self.pconfig = pconfig

    def get_headers_in_order(self) -> List[Tuple[int, str, Dict]]:
        """
        Gets the headers in the order they want to be displayed.
        Returns a list of triplets: (bucket_idx, key, header_info)
        """
        res = list()
        # Scan through self.headers_in_order and just bolt on the actual header info
        for bucket in sorted(self.headers_in_order):
            for bucket_idx, k in self.headers_in_order[bucket]:
                res.append((bucket_idx, k, self.headers[bucket_idx][k]))
        return res
