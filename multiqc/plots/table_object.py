"""
MultiQC datatable class, used by tables and violin plots
"""

import math

import logging
import re
from collections import defaultdict
from typing import List, Tuple, Dict, Optional, Union, Mapping, Callable

from pydantic import BaseModel, Field

from multiqc import config, report
from multiqc.plots.plotly.plot import PConfig

logger = logging.getLogger(__name__)


class TableConfig(PConfig):
    namespace: str = ""
    save_file: bool = False
    raw_data_fn: Optional[str] = None
    defaultsort: Optional[List[Dict[str, str]]] = None
    sort_rows: bool = Field(True, alias="sortRows")
    only_defined_headers: bool = True
    col1_header: str = "Sample Name"
    no_violin: bool = Field(False, alias="no_beeswarm")
    scale: Union[str, bool] = "GnBu"
    min: Optional[Union[int, float]] = None


class TableColumn(BaseModel):
    """
    Column model class. Holds configuration for a single column in a table.
    """

    rid: str
    title: str
    description: str
    namespace: str
    scale: Union[str, bool]
    hidden: bool
    color: str = Field(validation_alias="colour")
    placement: float = None
    max: Optional[float] = None
    dmax: Optional[float] = None
    min: Optional[float] = None
    dmin: Optional[float] = None
    ceiling: Optional[float] = None
    floor: Optional[float] = None
    minrange: Optional[float] = None
    shared_key: Optional[str]
    tt_decimals: Optional[int] = None
    suffix: Optional[str] = None
    cond_formatting_colours: List[Dict[str, str]] = []
    cond_formatting_rules: Dict[str, List[Dict[str, str]]] = {}
    bgcols: Dict[str, str] = {}
    bars_zero_centrepoint: bool = False


ValueT = Union[int, float, str, bool]


class DataTable(BaseModel):
    """
    Data table class. Prepares and holds data and configuration
    for either a table or a violin plot.
    """

    id: str
    raw_data: List[Dict[str, Dict[str, ValueT]]] = []
    formatted_data: List[Dict[str, Dict[str, str]]] = []
    headers_in_order: Dict[int, List[Tuple[int, str]]]
    headers: List[Dict[str, TableColumn]] = []
    pconfig: TableConfig

    @staticmethod
    def create(
        data: Union[List[Mapping[str, Mapping[str, Optional[ValueT]]]], Mapping[str, Mapping[str, Optional[ValueT]]]],
        pconfig: TableConfig,
        headers: Optional[Union[List[Dict[str, Dict]], Dict[str, Dict]]] = None,
    ) -> "DataTable":
        """Prepare data for use in a table or plot"""
        if headers is None:
            headers = []

        # Given one dataset - turn it into a list
        if not isinstance(data, list):
            data = [data]
        if not isinstance(headers, list):
            headers = [headers]

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

        raw_data = []
        formatted_data = []

        # Go through each table section
        for d_idx, dataset in enumerate(data):
            # Get the header keys
            try:
                keys = headers[d_idx].keys()
                assert len(keys) > 0
            except (IndexError, AttributeError, AssertionError):
                pconfig.only_defined_headers = False
                keys = list()

            # Add header keys from the data
            if not pconfig.only_defined_headers:
                # Get the keys from the data
                keys = list()
                for v_by_metric in dataset.values():
                    for k in v_by_metric.keys():
                        if k not in keys:
                            keys.append(k)

                # If we don't have a headers dict for this data set yet, create one
                try:
                    headers[d_idx]
                except IndexError:
                    headers.append(dict())
                else:
                    # Convert the existing headers into a dict (e.g. if parsed from a config)
                    od_tuples = [(key, headers[d_idx][key]) for key in headers[d_idx].keys()]
                    headers[d_idx] = dict(od_tuples)

                # Create empty header configs for each new data key
                for k in keys:
                    if k not in headers[d_idx]:
                        headers[d_idx][k] = {}

            # Ensure that keys are strings, not numeric
            keys = [str(k) for k in keys]
            for k in list(headers[d_idx].keys()):
                headers[d_idx][str(k)] = headers[d_idx].pop(k)

            # Ensure that all sample names are strings as well
            cdata = dict()
            for s_name, d in data[d_idx].items():
                cdata[str(s_name)] = d
            data[d_idx] = cdata

            # Ensure metric names are strings
            for s_name in data[d_idx].keys():
                for metric in list(data[d_idx][s_name].keys()):
                    data[d_idx][s_name][str(metric)] = data[d_idx][s_name].pop(metric)

            # Check that we have some data in each column
            empties = list()
            for k in keys:
                n = 0
                for v_by_metric in dataset.values():
                    if k in v_by_metric:
                        n += 1
                if n == 0:
                    empties.append(k)
            for k in empties:
                keys = [j for j in keys if j != k]
                logger.debug(
                    f"Table key '{k}' not found in data for '{pconfig.id}'. Skipping. Check for possible typos between data keys and header keys"
                )
                del headers[d_idx][k]

            raw_dataset: Dict[str, Dict[str, ValueT]] = defaultdict(dict)
            formatted_dataset: Dict[str, Dict[str, str]] = defaultdict(dict)
            for k in keys:
                # Unique id to avoid overwriting by other datasets
                if "rid" not in headers[d_idx][k]:
                    headers[d_idx][k]["rid"] = report.save_htmlid(re.sub(r"\W+", "_", k).strip().strip("_"))

                # Applying defaults presets for data keys if shared_key is set to base_count or read_count
                shared_key = headers[d_idx][k].get("shared_key", None)
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
                    if "modify" not in headers[d_idx][k]:
                        headers[d_idx][k]["modify"] = lambda x: x * multiplier
                    if "min" not in headers[d_idx][k] is None:
                        headers[d_idx][k]["min"] = 0
                    if "format" not in headers[d_idx][k] is None:
                        if multiplier == 1:
                            headers[d_idx][k]["format"] = "{:,d}"
                if "suffix" not in headers[d_idx][k] and shared_key_suffix is not None:
                    headers[d_idx][k]["suffix"] = " " + shared_key_suffix

                # Use defaults / data keys if headers not given
                headers[d_idx][k]["namespace"] = headers[d_idx][k].get("namespace", pconfig.namespace)
                headers[d_idx][k]["title"] = headers[d_idx][k].get("title", k)
                headers[d_idx][k]["description"] = headers[d_idx][k].get("description", headers[d_idx][k]["title"])
                headers[d_idx][k]["scale"] = headers[d_idx][k].get("scale", pconfig.scale)
                headers[d_idx][k]["format"] = headers[d_idx][k].get("format")
                headers[d_idx][k]["colour"] = headers[d_idx][k].get("colour", headers[d_idx][k].get("color"))
                headers[d_idx][k]["hidden"] = headers[d_idx][k].get("hidden", False)
                headers[d_idx][k]["max"] = headers[d_idx][k].get("max")
                headers[d_idx][k]["min"] = headers[d_idx][k].get("min", pconfig.min)
                headers[d_idx][k]["ceiling"] = headers[d_idx][k].get("ceiling")
                headers[d_idx][k]["floor"] = headers[d_idx][k].get("floor")
                headers[d_idx][k]["minrange"] = headers[d_idx][k].get("minrange", headers[d_idx][k].get("minRange"))
                headers[d_idx][k]["shared_key"] = headers[d_idx][k].get("shared_key")
                headers[d_idx][k]["modify"] = headers[d_idx][k].get("modify")
                headers[d_idx][k]["placement"] = float(headers[d_idx][k].get("placement", 1000))

                if headers[d_idx][k]["colour"] is None:
                    cidx = d_idx
                    while cidx >= len(SECTION_COLORS):
                        cidx -= len(SECTION_COLORS)
                    headers[d_idx][k]["colour"] = SECTION_COLORS[cidx]

                # Overwrite (2nd time) any given config with table-level user config
                # This is to override column-specific values set by modules
                if pconfig.id in config.custom_plot_config:
                    for cpc_k, cpc_v in config.custom_plot_config[pconfig.id].items():
                        headers[d_idx][k][cpc_k] = cpc_v

                # Overwrite "name" if set in user config
                # Key can be a column ID, a table ID, or a namespace in the general stats table.
                for key, val in config.table_columns_name.items():
                    key = key.lower()
                    # Case-insensitive check if the outer key is a table ID or a namespace.
                    if key in [pconfig.id.lower(), headers[d_idx][k]["namespace"].lower()] and isinstance(val, dict):
                        # Assume a dict of specific column IDs
                        for key2, new_title in val.items():
                            key2 = key2.lower()
                            if key2 in [k.lower(), headers[d_idx][k]["title"].lower()]:
                                headers[d_idx][k]["title"] = new_title

                    # Case-insensitive check if the outer key is a column ID
                    elif key in [k.lower(), headers[d_idx][k]["title"].lower()] and isinstance(val, str):
                        headers[d_idx][k]["title"] = val

                # Overwrite "hidden" if set in user config
                # Key can be a column ID, a table ID, or a namespace in the general stats table.
                for key, val in config.table_columns_visible.items():
                    key = key.lower()
                    # Case-insensitive check if the outer key is a table ID or a namespace.
                    if key in [pconfig.id.lower(), headers[d_idx][k]["namespace"].lower()]:
                        # First - if config value is a bool, set all module columns to that value
                        if isinstance(val, bool):
                            # Config has True = visible, False = Hidden. Here we're setting "hidden" which is inverse
                            headers[d_idx][k]["hidden"] = not val

                        # Not a bool, assume a dict of specific column IDs
                        elif isinstance(val, dict):
                            for key2, visible in val.items():
                                key2 = key2.lower()
                                if key2 in [k.lower(), headers[d_idx][k]["title"].lower()] and isinstance(
                                    visible, bool
                                ):
                                    # Config has True = visible, False = Hidden. Here we're setting "hidden" which is inverse
                                    headers[d_idx][k]["hidden"] = not visible

                    # Case-insensitive check if the outer key is a column ID
                    elif key in [k.lower(), headers[d_idx][k]["title"].lower()] and isinstance(val, bool):
                        # Config has True = visible, False = Hidden. Here we're setting "hidden" which is inverse
                        headers[d_idx][k]["hidden"] = not val

                # Also overwrite placement if set in config
                try:
                    headers[d_idx][k]["placement"] = float(
                        config.table_columns_placement[headers[d_idx][k]["namespace"]][k]
                    )
                except (KeyError, ValueError):
                    try:
                        headers[d_idx][k]["placement"] = float(config.table_columns_placement[pconfig.id][k])
                    except (KeyError, ValueError):
                        pass

                # Overwrite any header config if set in config
                for custom_k, custom_v in config.custom_table_header_config.get(pconfig.id, {}).get(k, {}).items():
                    headers[d_idx][k][custom_k] = custom_v

                # Filter keys and apply "modify" and "format" to values. Builds a copy of a dataset
                for s_name, v_by_metric in data[d_idx].items():
                    if k in v_by_metric:
                        val = v_by_metric[k]
                        if val is None or str(val).strip() == "":
                            continue

                        # Try parse as a number
                        if str(val).isdigit():
                            val = int(val)
                        else:
                            try:
                                val = float(val)
                            except ValueError:
                                pass

                        # Apply modify
                        if "modify" in headers[d_idx][k] and callable(headers[d_idx][k]["modify"]):
                            # noinspection PyBroadException
                            try:
                                val = headers[d_idx][k]["modify"](val)
                            except Exception as e:  # User-provided modify function can raise any exception
                                logger.error(f"Error modifying table value '{k}': '{val}'. {e}")
                        raw_dataset[s_name][k] = val

                        # Now also calculate formatted values
                        valstr = str(val)
                        fmt: Union[None, str, Callable] = headers[d_idx][k].get("format")
                        if fmt is None:
                            if isinstance(val, float):
                                fmt = "{:,.1f}"
                            elif isinstance(val, int):
                                fmt = "{:,d}"
                        if fmt is not None:
                            if callable(fmt):
                                try:
                                    # noinspection PyCallingNonCallable
                                    valstr = fmt(val)
                                except Exception as e:
                                    logger.error(f"Error applying format to table value '{k}': '{val}'. {e}")
                            elif isinstance(val, (int, float)):
                                try:
                                    valstr = fmt.format(val)
                                except Exception as e:
                                    logger.error(
                                        f"Error applying format string '{fmt}' to table value '{k}': '{val}'. {e}. "
                                        f"Check if your format string is correct."
                                    )
                        formatted_dataset[s_name][k] = valstr

                # Work out max and min value if not given
                setdmax = False
                setdmin = False
                try:
                    headers[d_idx][k]["dmax"] = float(headers[d_idx][k]["max"])
                except Exception:
                    headers[d_idx][k]["dmax"] = 0
                    setdmax = True

                try:
                    headers[d_idx][k]["dmin"] = float(headers[d_idx][k]["min"])
                except Exception:
                    headers[d_idx][k]["dmin"] = 0
                    setdmin = True

                # Figure out the min / max if not supplied
                if setdmax or setdmin:
                    for s_name, v_by_metric in raw_dataset.items():
                        try:
                            val = float(v_by_metric[k])
                            if math.isfinite(val) and not math.isnan(val):
                                if setdmax:
                                    headers[d_idx][k]["dmax"] = max(headers[d_idx][k]["dmax"], val)
                                if setdmin:
                                    headers[d_idx][k]["dmin"] = min(headers[d_idx][k]["dmin"], val)
                        except (ValueError, TypeError):
                            val = v_by_metric[k]  # couldn't convert to float - keep as a string
                        except KeyError:
                            pass  # missing data - skip
                    # Limit auto-generated scales with floor, ceiling and minrange.
                    if headers[d_idx][k]["ceiling"] is not None and headers[d_idx][k]["max"] is None:
                        headers[d_idx][k]["dmax"] = min(headers[d_idx][k]["dmax"], float(headers[d_idx][k]["ceiling"]))
                    if headers[d_idx][k]["floor"] is not None and headers[d_idx][k]["min"] is None:
                        headers[d_idx][k]["dmin"] = max(headers[d_idx][k]["dmin"], float(headers[d_idx][k]["floor"]))
                    if headers[d_idx][k]["minrange"] is not None:
                        drange = headers[d_idx][k]["dmax"] - headers[d_idx][k]["dmin"]
                        if drange < float(headers[d_idx][k]["minrange"]):
                            headers[d_idx][k]["dmax"] = headers[d_idx][k]["dmin"] + float(headers[d_idx][k]["minrange"])

            raw_data.append(raw_dataset)
            formatted_data.append(formatted_dataset)

        # Collect settings for shared keys
        shared_keys = defaultdict(lambda: dict())
        for d_idx, hs in enumerate(headers):
            for k in hs.keys():
                sk = headers[d_idx][k]["shared_key"]
                if sk is not None:
                    shared_keys[sk]["dmax"] = max(
                        headers[d_idx][k]["dmax"], shared_keys[sk].get("dmax", headers[d_idx][k]["dmax"])
                    )
                    shared_keys[sk]["dmin"] = min(
                        headers[d_idx][k]["dmin"], shared_keys[sk].get("dmin", headers[d_idx][k]["dmin"])
                    )

        # Overwrite shared key settings and at the same time assign to buckets for sorting
        # So the final ordering is:
        #   placement > section > explicit_ordering
        # Of course, the user can shuffle these manually.
        headers_in_order = defaultdict(list)
        for d_idx, hs in enumerate(headers):
            for k in hs.keys():
                sk = headers[d_idx][k]["shared_key"]
                if sk is not None:
                    headers[d_idx][k]["dmax"] = shared_keys[sk]["dmax"]
                    headers[d_idx][k]["dmin"] = shared_keys[sk]["dmin"]

                headers_in_order[headers[d_idx][k]["placement"]].append((d_idx, k))

        # # Skip any data that is not used in the table
        # # Would be ignored for making the table anyway, but can affect whether a violin plot is used
        # for d_idx, dataset in enumerate(raw_data):
        #     for s_name in list(dataset.keys()):
        #         if not any(h in data[d_idx][s_name].keys() for h in headers[d_idx]):
        #             del raw_data[d_idx][s_name]

        # Remove callable headers that are not compatible with JSON
        headers = [
            {metric: {k: v for k, v in d.items() if k not in ["modify", "format"]} for metric, d in h.items()}
            for h in headers
        ]

        # Assign to class
        return DataTable(
            id=pconfig.id,
            raw_data=raw_data,
            formatted_data=formatted_data,
            headers_in_order=dict(headers_in_order),
            headers=headers,
            pconfig=pconfig,
        )

    def get_headers_in_order(self) -> List[Tuple[int, str, TableColumn]]:
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
