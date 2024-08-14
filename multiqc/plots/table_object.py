"""
MultiQC datatable class, used by tables and violin plots
"""

import math

import logging
import re
from collections import defaultdict
from typing import List, Tuple, Dict, Optional, Union, Callable, Sequence, Mapping

from pydantic import BaseModel, Field

from multiqc import config, report
from multiqc.plots.plotly.plot import PConfig
from multiqc.validation import ValidatedConfig

logger = logging.getLogger(__name__)


class TableConfig(PConfig):
    namespace: str = ""
    save_file: bool = False
    raw_data_fn: Optional[str] = None
    defaultsort: Optional[List[Dict[str, str]]] = None
    sortRows: bool = Field(True, deprecated="sort_rows")
    sort_rows: bool = True
    only_defined_headers: bool = True
    col1_header: str = "Sample Name"
    no_beeswarm: bool = Field(False, deprecated="no_violin")
    no_violin: bool = False
    scale: Union[str, bool] = "GnBu"
    min: Optional[Union[int, float]] = None


class TableColumn(ValidatedConfig):
    """
    Column model class. Holds configuration for a single column in a table.
    """

    id: Optional[str] = Field(None, deprecated="rid")
    rid: str
    title: str
    description: str
    namespace: str
    scale: Union[str, bool]
    hidden: bool
    colour: Optional[str] = Field(None, deprecated="color")
    color: Optional[str] = None
    placement: Optional[float] = None
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

DatasetT = Mapping[str, Union[Mapping[str, Optional[ValueT]], List[Tuple[str, Mapping[str, Optional[ValueT]]]]]]


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
        data: Union[DatasetT, Sequence[DatasetT]],
        pconfig: TableConfig,
        headers: Optional[Union[List[Dict[str, Dict]], Dict[str, Dict]]] = None,
    ) -> "DataTable":
        """Prepare data for use in a table or plot"""
        if headers is None:
            headers = []

        # Given one dataset - turn it into a list
        datasets__with_nulls = data if isinstance(data, list) else [data]
        list_of_headers = headers if isinstance(headers, list) else [headers]
        del data
        del headers

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

        raw_data: List[Dict[str, Dict[str, ValueT]]] = []
        formatted_data: List[Dict[str, Dict[str, str]]] = []

        # Go through each table section
        for d_idx, ds__with_nulls in enumerate(datasets__with_nulls):
            # Get the header keys
            if d_idx < len(list_of_headers):
                col_ids = list(list_of_headers[d_idx].keys())
                assert len(col_ids) > 0
            else:
                pconfig.only_defined_headers = False
                col_ids = list()

            # Add header keys from the data
            if not pconfig.only_defined_headers:
                # Get the keys from the data
                col_ids = list()
                for v_by_metric__with_nulls in ds__with_nulls.values():
                    for col_id in v_by_metric__with_nulls.keys():
                        if col_id not in col_ids:
                            col_ids.append(col_id)

                # If we don't have a headers dict for this data set yet, create one
                try:
                    list_of_headers[d_idx]
                except IndexError:
                    list_of_headers.append(dict())
                else:
                    # Convert the existing headers into a dict (e.g. if parsed from a config)
                    od_tuples = [(key, list_of_headers[d_idx][key]) for key in list_of_headers[d_idx].keys()]
                    list_of_headers[d_idx] = dict(od_tuples)

                # Create empty header configs for each new data key
                for col_id in col_ids:
                    if col_id not in list_of_headers[d_idx]:
                        list_of_headers[d_idx][col_id] = {}

            # Check that we have some data in each column
            empties = list()
            for col_id in col_ids:
                n = 0
                for v_by_metric__with_nulls in ds__with_nulls.values():
                    if isinstance(v_by_metric__with_nulls, dict):
                        if col_id in v_by_metric__with_nulls:
                            n += 1
                    elif isinstance(v_by_metric__with_nulls, list):
                        for _, v_by_metric in v_by_metric__with_nulls:
                            if col_id in v_by_metric:
                                n += 1
                if n == 0:
                    empties.append(col_id)
            for col_id in empties:
                col_ids = [j for j in col_ids if j != col_id]
                logger.debug(
                    f"Table key '{col_id}' not found in data for '{pconfig.id}'. Skipping. Check for possible typos between data keys and header keys"
                )
                del list_of_headers[d_idx][col_id]

            raw_dataset: Dict[str, Dict[str, ValueT]] = defaultdict(dict)
            formatted_dataset: Dict[str, Dict[str, str]] = defaultdict(dict)
            for col_id in col_ids:
                # Unique id to avoid overwriting by other datasets
                unclean_rid = list_of_headers[d_idx][col_id].get("rid", col_id)
                rid = re.sub(r"\W+", "_", unclean_rid).strip().strip("_")
                list_of_headers[d_idx][col_id]["rid"] = report.save_htmlid(report.clean_htmlid(rid), skiplint=True)

                # Applying defaults presets for data keys if shared_key is set to base_count or read_count
                shared_key = list_of_headers[d_idx][col_id].get("shared_key", None)
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
                    if "modify" not in list_of_headers[d_idx][col_id]:
                        list_of_headers[d_idx][col_id]["modify"] = lambda x: x * multiplier
                    if "min" not in list_of_headers[d_idx][col_id] is None:
                        list_of_headers[d_idx][col_id]["min"] = 0
                    if "format" not in list_of_headers[d_idx][col_id] is None:
                        if multiplier == 1:
                            list_of_headers[d_idx][col_id]["format"] = "{:,d}"
                if "suffix" not in list_of_headers[d_idx][col_id] and shared_key_suffix is not None:
                    list_of_headers[d_idx][col_id]["suffix"] = " " + shared_key_suffix

                # Use defaults / data keys if headers not given
                list_of_headers[d_idx][col_id]["namespace"] = list_of_headers[d_idx][col_id].get(
                    "namespace", pconfig.namespace
                )
                list_of_headers[d_idx][col_id]["title"] = list_of_headers[d_idx][col_id].get("title", col_id)
                list_of_headers[d_idx][col_id]["description"] = list_of_headers[d_idx][col_id].get(
                    "description", list_of_headers[d_idx][col_id]["title"]
                )
                list_of_headers[d_idx][col_id]["scale"] = list_of_headers[d_idx][col_id].get("scale", pconfig.scale)
                list_of_headers[d_idx][col_id]["format"] = list_of_headers[d_idx][col_id].get("format")
                list_of_headers[d_idx][col_id]["color"] = list_of_headers[d_idx][col_id].get(
                    "colour", list_of_headers[d_idx][col_id].get("color")
                )
                list_of_headers[d_idx][col_id]["hidden"] = list_of_headers[d_idx][col_id].get("hidden", False)
                list_of_headers[d_idx][col_id]["max"] = list_of_headers[d_idx][col_id].get("max")
                list_of_headers[d_idx][col_id]["min"] = list_of_headers[d_idx][col_id].get("min", pconfig.min)
                list_of_headers[d_idx][col_id]["ceiling"] = list_of_headers[d_idx][col_id].get("ceiling")
                list_of_headers[d_idx][col_id]["floor"] = list_of_headers[d_idx][col_id].get("floor")
                list_of_headers[d_idx][col_id]["minrange"] = list_of_headers[d_idx][col_id].get(
                    "minrange", list_of_headers[d_idx][col_id].get("minRange")
                )
                list_of_headers[d_idx][col_id]["shared_key"] = list_of_headers[d_idx][col_id].get("shared_key")
                list_of_headers[d_idx][col_id]["modify"] = list_of_headers[d_idx][col_id].get("modify")
                list_of_headers[d_idx][col_id]["placement"] = float(
                    list_of_headers[d_idx][col_id].get("placement", 1000)
                )

                if list_of_headers[d_idx][col_id]["color"] is None:
                    cidx = d_idx
                    while cidx >= len(SECTION_COLORS):
                        cidx -= len(SECTION_COLORS)
                    list_of_headers[d_idx][col_id]["color"] = SECTION_COLORS[cidx]

                # Overwrite (2nd time) any given config with table-level user config
                # This is to override column-specific values set by modules
                if pconfig.id in config.custom_plot_config:
                    for cpc_k, cpc_v in config.custom_plot_config[pconfig.id].items():
                        if cpc_k in TableColumn.model_fields.keys():
                            list_of_headers[d_idx][col_id][cpc_k] = cpc_v

                # Overwrite "name" if set in user config
                # Key can be a column ID, a table ID, or a namespace in the general stats table.
                for item_id, new_title_val in config.table_columns_name.items():
                    item_id = item_id.lower()
                    # Case-insensitive check if the outer key is a table ID or a namespace.
                    if item_id in [
                        pconfig.id.lower(),
                        list_of_headers[d_idx][col_id]["namespace"].lower(),
                    ] and isinstance(new_title_val, dict):
                        # Assume a dict of specific column IDs
                        for item_id2, new_title in new_title_val.items():
                            item_id2 = item_id2.lower()
                            if item_id2 in [col_id.lower(), list_of_headers[d_idx][col_id]["title"].lower()]:
                                list_of_headers[d_idx][col_id]["title"] = new_title

                    # Case-insensitive check if the outer key is a column ID
                    elif item_id in [col_id.lower(), list_of_headers[d_idx][col_id]["title"].lower()] and isinstance(
                        new_title_val, str
                    ):
                        list_of_headers[d_idx][col_id]["title"] = new_title_val

                # Overwrite "hidden" if set in user config
                # Key can be a column ID, a table ID, or a namespace in the general stats table.
                for item_id, visibility in config.table_columns_visible.items():
                    item_id = item_id.lower()
                    # Case-insensitive check if the outer key is a table ID or a namespace.
                    if item_id in [pconfig.id.lower(), list_of_headers[d_idx][col_id]["namespace"].lower()]:
                        # First - if config value is a bool, set all module columns to that value
                        if isinstance(visibility, bool):
                            # Config has True = visible, False = Hidden. Here we're setting "hidden" which is inverse
                            list_of_headers[d_idx][col_id]["hidden"] = not visibility

                        # Not a bool, assume a dict of specific column IDs
                        elif isinstance(visibility, dict):
                            for item_id2, visible in visibility.items():
                                item_id2 = item_id2.lower()
                                if item_id2 in [
                                    col_id.lower(),
                                    list_of_headers[d_idx][col_id]["title"].lower(),
                                ] and isinstance(visible, bool):
                                    # Config has True = visible, False = Hidden. Here we're setting "hidden" which is inverse
                                    list_of_headers[d_idx][col_id]["hidden"] = not visible

                    # Case-insensitive check if the outer key is a column ID
                    elif visibility in [col_id.lower(), list_of_headers[d_idx][col_id]["title"].lower()] and isinstance(
                        visibility, bool
                    ):
                        # Config has True = visible, False = Hidden. Here we're setting "hidden" which is inverse
                        list_of_headers[d_idx][col_id]["hidden"] = not visibility

                # Also overwrite placement if set in config
                try:
                    list_of_headers[d_idx][col_id]["placement"] = float(
                        config.table_columns_placement[list_of_headers[d_idx][col_id]["namespace"]][col_id]
                    )
                except (KeyError, ValueError):
                    try:
                        list_of_headers[d_idx][col_id]["placement"] = float(
                            config.table_columns_placement[pconfig.id][col_id]
                        )
                    except (KeyError, ValueError):
                        pass

                # Overwrite any header config if set in config
                for custom_k, custom_v in config.custom_table_header_config.get(pconfig.id, {}).get(col_id, {}).items():
                    list_of_headers[d_idx][col_id][custom_k] = custom_v

                # Filter keys and apply "modify" and "format" to values. Builds a copy of a dataset
                for s_name, v_by_metric__with_nulls in datasets__with_nulls[d_idx].items():
                    if col_id in v_by_metric__with_nulls:
                        val__with_nulls = v_by_metric__with_nulls[col_id]
                        if val__with_nulls is None or str(val__with_nulls).strip() == "":
                            continue
                        val: ValueT = val__with_nulls

                        # Try parse as a number
                        if str(val).isdigit():
                            val = int(val)
                        else:
                            try:
                                val = float(val)
                            except ValueError:
                                pass

                        # Apply modify
                        if "modify" in list_of_headers[d_idx][col_id] and callable(
                            list_of_headers[d_idx][col_id]["modify"]
                        ):
                            # noinspection PyBroadException
                            try:
                                val = list_of_headers[d_idx][col_id]["modify"](val)
                            except Exception as e:  # User-provided modify function can raise any exception
                                logger.error(f"Error modifying table value '{col_id}': '{val}'. {e}")
                        raw_dataset[s_name][col_id] = val

                        # Now also calculate formatted values
                        valstr = str(val)
                        fmt: Union[None, str, Callable] = list_of_headers[d_idx][col_id].get("format")
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
                                    logger.error(f"Error applying format to table value '{col_id}': '{val}'. {e}")
                            elif isinstance(val, (int, float)):
                                try:
                                    valstr = fmt.format(val)
                                except Exception as e:
                                    logger.error(
                                        f"Error applying format string '{fmt}' to table value '{col_id}': '{val}'. {e}. "
                                        f"Check if your format string is correct."
                                    )
                        formatted_dataset[s_name][col_id] = valstr

                # Work out max and min value if not given
                set_dmax = False
                set_dmin = False
                try:
                    list_of_headers[d_idx][col_id]["dmax"] = float(list_of_headers[d_idx][col_id]["max"])
                except Exception:
                    list_of_headers[d_idx][col_id]["dmax"] = 0
                    set_dmax = True

                try:
                    list_of_headers[d_idx][col_id]["dmin"] = float(list_of_headers[d_idx][col_id]["min"])
                except Exception:
                    list_of_headers[d_idx][col_id]["dmin"] = 0
                    set_dmin = True

                # Figure out the min / max if not supplied
                if set_dmax or set_dmin:
                    for s_name, v_by_metric in raw_dataset.items():
                        try:
                            val = float(v_by_metric[col_id])
                            if math.isfinite(val) and not math.isnan(val):
                                if set_dmax:
                                    list_of_headers[d_idx][col_id]["dmax"] = max(
                                        list_of_headers[d_idx][col_id]["dmax"], val
                                    )
                                if set_dmin:
                                    list_of_headers[d_idx][col_id]["dmin"] = min(
                                        list_of_headers[d_idx][col_id]["dmin"], val
                                    )
                        except (ValueError, TypeError):
                            val = v_by_metric[col_id]  # couldn't convert to float - keep as a string
                        except KeyError:
                            pass  # missing data - skip
                    # Limit auto-generated scales with floor, ceiling and minrange.
                    if (
                        list_of_headers[d_idx][col_id]["ceiling"] is not None
                        and list_of_headers[d_idx][col_id]["max"] is None
                    ):
                        list_of_headers[d_idx][col_id]["dmax"] = min(
                            list_of_headers[d_idx][col_id]["dmax"], float(list_of_headers[d_idx][col_id]["ceiling"])
                        )
                    if (
                        list_of_headers[d_idx][col_id]["floor"] is not None
                        and list_of_headers[d_idx][col_id]["min"] is None
                    ):
                        list_of_headers[d_idx][col_id]["dmin"] = max(
                            list_of_headers[d_idx][col_id]["dmin"], float(list_of_headers[d_idx][col_id]["floor"])
                        )
                    if list_of_headers[d_idx][col_id]["minrange"] is not None:
                        drange = list_of_headers[d_idx][col_id]["dmax"] - list_of_headers[d_idx][col_id]["dmin"]
                        if drange < float(list_of_headers[d_idx][col_id]["minrange"]):
                            list_of_headers[d_idx][col_id]["dmax"] = list_of_headers[d_idx][col_id]["dmin"] + float(
                                list_of_headers[d_idx][col_id]["minrange"]
                            )

            raw_data.append(raw_dataset)
            formatted_data.append(formatted_dataset)

        # Collect settings for shared keys
        shared_keys: Dict[str, Dict[str, Union[int, float]]] = defaultdict(lambda: dict())
        for d_idx, hs in enumerate(list_of_headers):
            for col_id in hs.keys():
                sk: str = list_of_headers[d_idx][col_id]["shared_key"]
                if sk is not None:
                    shared_keys[sk]["dmax"] = max(
                        list_of_headers[d_idx][col_id]["dmax"],
                        shared_keys[sk].get("dmax", list_of_headers[d_idx][col_id]["dmax"]),
                    )
                    shared_keys[sk]["dmin"] = min(
                        list_of_headers[d_idx][col_id]["dmin"],
                        shared_keys[sk].get("dmin", list_of_headers[d_idx][col_id]["dmin"]),
                    )

        # Overwrite shared key settings and at the same time assign to buckets for sorting
        # So the final ordering is:
        #   placement > section > explicit_ordering
        # Of course, the user can shuffle these manually.
        headers_in_order = defaultdict(list)
        for d_idx, hs in enumerate(list_of_headers):
            for col_id in hs.keys():
                sk = list_of_headers[d_idx][col_id]["shared_key"]
                if sk is not None:
                    list_of_headers[d_idx][col_id]["dmax"] = shared_keys[sk]["dmax"]
                    list_of_headers[d_idx][col_id]["dmin"] = shared_keys[sk]["dmin"]

                headers_in_order[list_of_headers[d_idx][col_id]["placement"]].append((d_idx, col_id))

        # # Skip any data that is not used in the table
        # # Would be ignored for making the table anyway, but can affect whether a violin plot is used
        # for d_idx, dataset in enumerate(raw_data):
        #     for s_name in list(dataset.keys()):
        #         if not any(h in data[d_idx][s_name].keys() for h in headers[d_idx]):
        #             del raw_data[d_idx][s_name]

        # Remove callable headers that are not compatible with JSON
        list_of_headers = [
            {metric: {k: v for k, v in d.items() if k not in ["modify", "format"]} for metric, d in h.items()}
            for h in list_of_headers
        ]

        # Assign to class
        return DataTable(
            id=pconfig.id,
            raw_data=raw_data,
            formatted_data=formatted_data,
            headers_in_order=dict(headers_in_order),
            headers=list_of_headers,
            pconfig=pconfig,
        )

    def get_headers_in_order(self) -> List[Tuple[int, str, TableColumn]]:
        """
        Gets the headers in the order they want to be displayed.
        Returns a list of triplets: (bucket_idx, key, header_info)
        """
        res: List[Tuple[int, str, TableColumn]] = list()
        # Scan through self.headers_in_order and just bolt on the actual header info
        for bucket in sorted(self.headers_in_order):
            for bucket_idx, k in self.headers_in_order[bucket]:
                res.append((bucket_idx, k, self.headers[bucket_idx][k]))
        return res
