"""
MultiQC datatable class, used by tables and violin plots
"""

import logging
import math
import re
from collections import defaultdict
from typing import Callable, Dict, List, Mapping, NewType, Optional, Sequence, Set, Tuple, Union

from pydantic import BaseModel, Field

from multiqc import config, report
from multiqc.plots.plotly.plot import PConfig
from multiqc.types import AnchorT, ColumnKeyT, SampleGroupT, SampleNameT
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


ColumnAnchorT = NewType("ColumnAnchorT", str)  # Unique within a table


class ColumnMeta(ValidatedConfig):
    """
    Column model class. Holds configuration for a single column in a table.
    """

    rid: ColumnAnchorT  # namespace + short_rid = ID unique within a table
    title: str
    description: str
    scale: Union[str, bool]
    hidden: bool
    placement: float
    namespace: str = ""
    colour: Optional[str] = Field(None, deprecated="color")
    color: Optional[str] = None
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
    modify: Optional[Callable] = None
    format: Optional[Union[str, Callable]] = None

    @staticmethod
    def create(
        header_d: Dict[str, Union[str, int, float, None, Callable]],
        col_key: ColumnKeyT,  # to initialize rid
        sec_idx: int,  # to initialize the colour
        pconfig: TableConfig,  # plot config dictionary
        table_anchor: AnchorT,
    ) -> "ColumnMeta":
        ns = header_d.get("namespace", pconfig.namespace) or ""
        assert isinstance(ns, str)
        if ns:
            header_d["namespace"] = ns

        unclean_rid = header_d.get("rid") or col_key
        legacy_short_rid = re.sub(r"\W+", "_", str(unclean_rid)).strip().strip("_")

        # Unique id to avoid overwriting by other datasets
        rid = legacy_short_rid
        if ns:
            ns = re.sub(r"\W+", "_", str(ns)).strip().strip("_").lower()
            rid = f"{ns}-{rid}"
        header_d["rid"] = ColumnAnchorT(report.save_htmlid(rid, scope=table_anchor))

        # Applying defaults presets for data keys if shared_key is set to base_count or read_count
        shared_key = header_d.get("shared_key", None)
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
            if "modify" not in header_d:
                header_d["modify"] = lambda x: x * multiplier
            if "min" not in header_d is None:
                header_d["min"] = 0
            if "format" not in header_d is None:
                if multiplier == 1:
                    header_d["format"] = "{:,d}"
        if "suffix" not in header_d and shared_key_suffix is not None:
            header_d["suffix"] = " " + shared_key_suffix

        # Use defaults / data keys if headers not given
        header_d["title"] = header_d.get("title", col_key)
        header_d["description"] = header_d.get("description", header_d["title"])
        header_d["scale"] = header_d.get("scale", pconfig.scale)
        header_d["format"] = header_d.get("format")
        header_d["color"] = header_d.get("colour", header_d.get("color"))
        header_d["hidden"] = header_d.get("hidden", False)
        header_d["max"] = header_d.get("max")
        header_d["min"] = header_d.get("min", pconfig.min)
        header_d["ceiling"] = header_d.get("ceiling")
        header_d["floor"] = header_d.get("floor")
        header_d["minrange"] = header_d.get("minrange", header_d.get("minRange"))
        header_d["shared_key"] = header_d.get("shared_key")
        modify = header_d.get("modify")
        if modify is not None and callable(modify):
            header_d["modify"] = modify
        placement = header_d.get("placement")
        if placement is not None and isinstance(placement, (str, float, int)):
            header_d["placement"] = float(placement)
        else:
            header_d["placement"] = 1000  # default value

        if header_d["color"] is None:
            cidx = sec_idx
            while cidx >= len(SECTION_COLORS):
                cidx -= len(SECTION_COLORS)
            header_d["color"] = SECTION_COLORS[cidx]

        # Overwrite (2nd time) any given config with table-level user config
        # This is to override column-specific values set by modules
        if pconfig.id in config.custom_plot_config:
            for cpc_k, cpc_v in config.custom_plot_config[pconfig.id].items():
                if cpc_k in ColumnMeta.model_fields.keys():
                    header_d[cpc_k] = cpc_v

        col: ColumnMeta = ColumnMeta(**header_d)

        def _ns_match(item_id: str) -> bool:
            return item_id.lower() in [
                str(s).lower()
                for s in [
                    pconfig.id,
                    pconfig.anchor,
                    col.namespace,
                ]
                if s is not None
            ]

        def _col_match(item_id: str) -> bool:
            return item_id.lower() in [
                str(s).lower()
                for s in [
                    col.rid,
                    legacy_short_rid,
                    col_key,
                    col.title,
                ]
                if s is not None
            ]

        # Overwrite "name" if set in user config
        # Key can be a column ID, a table ID, or a namespace in the general stats table.
        for item_id, new_title_val in config.table_columns_name.items():
            # Case-insensitive check if the outer key is a table ID or a namespace.
            if _ns_match(item_id) and isinstance(new_title_val, dict):
                # Assume a dict of specific column IDs
                for item_id2, new_title in new_title_val.items():
                    item_id2 = item_id2.lower()
                    if _col_match(item_id2):
                        col.title = new_title

            # Case-insensitive check if the outer key is a column ID
            elif _col_match(item_id) and isinstance(new_title_val, str):
                col.title = new_title_val

        # Overwrite "hidden" if set in user config
        # Key can be a column ID, a table ID, or a namespace in the general stats table.
        for item_id, visibility in config.table_columns_visible.items():
            # Case-insensitive check if the outer key is a table ID or a namespace.
            if _ns_match(item_id):
                # First - if config value is a bool, set all module columns to that value
                if isinstance(visibility, bool):
                    # Config has True = visible, False = Hidden. Here we're setting "hidden" which is inverse
                    col.hidden = not visibility

                # Not a bool, assume a dict of specific column IDs
                elif isinstance(visibility, dict):
                    for item_id2, visible in visibility.items():
                        item_id2 = item_id2.lower()
                        if _col_match(item_id2) and isinstance(visible, bool):
                            # Config has True = visible, False = Hidden. Here we're setting "hidden" which is inverse
                            col.hidden = not visible

            # Case-insensitive check if the outer key is a column ID
            elif _col_match(item_id) and isinstance(visibility, bool):
                # Config has True = visible, False = Hidden. Here we're setting "hidden" which is inverse
                col.hidden = not visibility

        # Also overwrite placement if set in config
        try:
            col.placement = float(config.table_columns_placement[col.namespace][str(col_key)])
        except (KeyError, ValueError):
            try:
                col.placement = float(config.table_columns_placement[pconfig.id][str(col_key)])
            except (KeyError, ValueError):
                pass

        # Overwrite any header config if set in config
        for custom_k, custom_v in config.custom_table_header_config.get(pconfig.id, {}).get(col_key, {}).items():
            setattr(col, custom_k, custom_v)

        return col


ValueT = Union[int, float, str, bool]


class InputRowT(BaseModel):
    """
    Row class. Holds configuration for a single row in a table (can be multiple for one sample)
    """

    sample: SampleNameT
    data: Dict[ColumnKeyT, Optional[ValueT]] = dict()


ColumnKey = Union[str, ColumnKeyT]
SampleGroup = Union[str, SampleGroupT]
InputGroupT = Union[Mapping[ColumnKey, Optional[ValueT]], InputRowT, Sequence[InputRowT]]
InputSectionT = Mapping[SampleGroup, InputGroupT]
InputHeaderT = Mapping[ColumnKey, Mapping[str, Union[str, int, float, None, Callable]]]


class Row(BaseModel):
    """
    Processed row class. Holds configuration for a single row in a table (can be multiple for one sample).
    Contains raw, optionally modified, non-null values, and corresponding formatted values to display.
    """

    sample: SampleNameT
    # rows with original, unformatted values coming from modules:
    raw_data: Dict[ColumnKeyT, ValueT] = dict()
    # formatted rows (i.e. values are HTML strings to display in the table):
    formatted_data: Dict[ColumnKeyT, str] = dict()


class TableSection(BaseModel):
    """
    Table section class. Holds configuration for a single section in a table.
    """

    column_by_key: Dict[ColumnKeyT, ColumnMeta]
    rows_by_sgroup: Dict[SampleGroupT, List[Row]] = defaultdict(list)


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


col_anchors_by_table: Dict[AnchorT, Set[ColumnAnchorT]] = defaultdict(set)


class DataTable(BaseModel):
    """
    Data table class. Prepares and holds data and configuration
    for either a table or a violin plot.
    """

    id: str
    anchor: AnchorT
    pconfig: TableConfig

    sections: List[TableSection]
    headers_in_order: Dict[float, List[Tuple[int, ColumnKeyT]]]

    @staticmethod
    def create(
        data: Union[InputSectionT, Sequence[InputSectionT]],
        table_id: str,
        table_anchor: AnchorT,
        pconfig: TableConfig,
        headers: Optional[Union[List[InputHeaderT], InputHeaderT]] = None,
    ) -> "DataTable":
        """Prepare data for use in a table or plot"""
        if headers is None:
            headers = []

        # Given one dataset - turn it into a list
        input_sections__with_nulls = data if isinstance(data, list) else [data]
        list_of_headers = headers if isinstance(headers, list) else [headers]
        del data
        del headers

        # Each section to have a list of groups (even if there is just one element in a group)
        input_section: InputSectionT
        input_group: InputGroupT
        unified_sections__with_nulls: List[Dict[SampleGroupT, List[InputRowT]]] = []
        for input_section in input_sections__with_nulls:
            rows_by_group: Dict[SampleGroupT, List[InputRowT]] = {}
            for g_name, input_group in input_section.items():
                g_name = SampleGroupT(str(g_name))  # Make sure sample names are strings
                if isinstance(input_group, dict):  # just one row, defined as a mapping from metric to value
                    # Remove non-scalar values for table cells
                    input_group = {k: v for k, v in input_group.items() if isinstance(v, (int, float, str, bool))}
                    rows_by_group[g_name] = [InputRowT(sample=g_name, data=input_group)]
                elif isinstance(input_group, list):  # multiple rows, each defined as a mapping from metric to value
                    rows_by_group[g_name] = input_group
                else:
                    assert isinstance(input_group, InputRowT)
                    rows_by_group[g_name] = [input_group]
            unified_sections__with_nulls.append(rows_by_group)

        del input_sections__with_nulls

        # Go through each table section and create a list of Section objects
        sections: List[TableSection] = []
        for sec_idx, rows_by_sname__with_nulls in enumerate(unified_sections__with_nulls):
            header_by_key: InputHeaderT = list_of_headers[sec_idx] if sec_idx < len(list_of_headers) else dict()
            if not header_by_key:
                pconfig.only_defined_headers = False

            column_by_key: Dict[ColumnKeyT, ColumnMeta] = dict()
            header_by_key_copy = _get_or_create_headers(rows_by_sname__with_nulls, header_by_key, pconfig)
            for col_key, header_d in header_by_key_copy.items():
                column_by_key[col_key] = ColumnMeta.create(
                    header_d=header_d, col_key=col_key, sec_idx=sec_idx, pconfig=pconfig, table_anchor=table_anchor
                )

            # Filter out null values and columns that are not present in column_by_key,
            # and apply "modify" and "format" to values. Will generate non-null data and str data.
            section = TableSection(column_by_key=column_by_key)
            for g_name, group_rows__with_nulls in rows_by_sname__with_nulls.items():
                for input_row in group_rows__with_nulls:
                    row = Row(sample=input_row.sample)
                    for col_key, optional_val in input_row.data.items():
                        if col_key not in column_by_key:  # missing in provided headers
                            continue
                        if optional_val is None or str(optional_val).strip() == "":  # empty
                            continue
                        val, valstr = _process_and_format_value(optional_val, column_by_key[col_key])
                        row.raw_data[col_key] = val
                        row.formatted_data[col_key] = valstr
                    if row.raw_data:
                        section.rows_by_sgroup[g_name].append(row)

            # Remove empty groups:
            section.rows_by_sgroup = {sname: rows for sname, rows in section.rows_by_sgroup.items() if rows}

            # Work out max and min value if not given:
            for col_key, column in column_by_key.items():
                _determine_dmin_and_dmax(column, col_key, section.rows_by_sgroup)

            sections.append(section)

        del unified_sections__with_nulls
        del list_of_headers

        shared_keys: Dict[str, Dict[str, Union[int, float]]] = _collect_shared_keys(sections)

        # Overwrite shared key settings and at the same time assign to buckets for sorting
        # So the final ordering is:
        #   placement > section > explicit_ordering
        # Of course, the user can shuffle these manually.
        headers_in_order: Dict[float, List[Tuple[int, ColumnKeyT]]] = defaultdict(list)
        for sec_idx, section in enumerate(sections):
            for col_key, column in section.column_by_key.items():
                if column.shared_key is not None:
                    column.dmax = shared_keys[column.shared_key]["dmax"]
                    column.dmin = shared_keys[column.shared_key]["dmin"]

                headers_in_order[column.placement].append((sec_idx, col_key))

        # # Skip any data that is not used in the table
        # # Would be ignored for making the table anyway, but can affect whether a violin plot is used
        # for d_idx, dataset in enumerate(raw_data):
        #     for s_name in list(dataset.keys()):
        #         if not any(h in data[d_idx][s_name].keys() for h in headers[d_idx]):
        #             del raw_data[d_idx][s_name]

        # Remove callable headers that are not compatible with JSON
        for section in sections:
            for column in section.column_by_key.values():
                del column.modify
                del column.format

        # Assign to class
        return DataTable(
            id=table_id,
            anchor=table_anchor,
            sections=sections,
            headers_in_order=headers_in_order,
            pconfig=pconfig,
        )

    def get_headers_in_order(self) -> List[Tuple[int, ColumnKeyT, ColumnMeta]]:
        """
        Gets the headers in the order they want to be displayed.
        Returns a list of triplets: (bucket_idx, key, header_info)
        """
        res: List[Tuple[int, ColumnKeyT, ColumnMeta]] = list()
        # Scan through self.headers_in_order and just bolt on the actual header info
        placement: float
        for placement in sorted(self.headers_in_order.keys()):
            for section_idx, col_key in self.headers_in_order[placement]:
                res.append((section_idx, col_key, self.sections[section_idx].column_by_key[col_key]))
        return res


def _get_or_create_headers(
    rows_by_sample: Dict[SampleGroupT, List[InputRowT]],
    header_by_key: InputHeaderT,
    pconfig,
) -> Dict[ColumnKeyT, Dict[str, Union[str, int, float, None, Callable]]]:
    """
    Process and populate headers, if missing or incomplete.
    """
    # Make a copy to keep the input immutable.
    header_by_key_copy = {ColumnKeyT(k): dict(h) for k, h in header_by_key.items()}
    if not pconfig.only_defined_headers:
        # Get additional header keys from the data
        col_ids: List[ColumnKeyT] = list(header_by_key_copy.keys())
        # Get the keys from the data
        for sname, rows in rows_by_sample.items():
            for row in rows:
                for col_id in row.data.keys():
                    if col_id not in col_ids:
                        col_ids.append(col_id)

        # Create empty header configs for each new data key
        for col_id in col_ids:
            if col_id not in header_by_key:
                header_by_key_copy[col_id] = dict()

    # Check that we have some data in each column
    empties = list()
    for col_id in header_by_key_copy.keys():
        n = 0
        for sname, rows in rows_by_sample.items():
            for row in rows:
                if col_id in row.data.keys():
                    n += 1

        if n == 0:
            empties.append(col_id)

    # Remove empty columns
    for empty_col_id in empties:
        logger.debug(
            f"Table key '{empty_col_id}' not found in data for '{pconfig.id}'. Skipping. Check for possible typos between data keys and header keys"
        )
        del header_by_key_copy[empty_col_id]

    return header_by_key_copy


def _process_and_format_value(val: ValueT, column: ColumnMeta) -> Tuple[ValueT, str]:
    """
    Takes row value, applies "modify" and "format" functions, and returns a tuple:
    the modified value and its formatted string.
    """
    # Try parse as a number
    if str(val).isdigit():
        val = int(val)
    else:
        try:
            val = float(val)
        except ValueError:
            pass

    # Apply modify
    if column.modify:
        # noinspection PyBroadException
        try:
            val = column.modify(val)
        except Exception as e:  # User-provided modify function can raise any exception
            logger.error(f"Error modifying table value '{column.id}': '{val}'. {e}")
    # section.raw_data[s_name][col_key] = val

    # Now also calculate formatted values
    valstr = str(val)
    fmt: Union[None, str, Callable] = column.format
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
                logger.error(f"Error applying format to table value '{column.rid}': '{val}'. {e}")
        elif isinstance(val, (int, float)):
            try:
                valstr = fmt.format(val)
            except Exception as e:
                logger.error(
                    f"Error applying format string '{fmt}' to table value '{column.rid}': '{val}'. {e}. "
                    f"Check if your format string is correct."
                )
    # section.formatted_data[s_name][col_key] = valstr
    return val, valstr


def _determine_dmin_and_dmax(
    column: ColumnMeta,
    col_key: ColumnKeyT,
    rows_by_sample: Dict[SampleGroupT, List[Row]],
) -> None:
    """
    Work out max and min value in a column if not given, to support color scale.
    """

    set_dmax = False
    if column.max is not None:
        column.dmax = float(column.max)
    else:
        column.dmax = 0
        set_dmax = True

    set_dmin = False
    if column.min is not None:
        column.dmin = float(column.min)
    else:
        column.dmin = 0
        set_dmin = True

    # Figure out the min / max if not supplied
    if set_dmax or set_dmin:
        for s_name, rows in rows_by_sample.items():
            v_by_col = rows[0].raw_data
            try:
                val = v_by_col[col_key]
                if isinstance(val, float) or isinstance(val, int) and math.isfinite(val) and not math.isnan(val):
                    if set_dmax:
                        column.dmax = max(column.dmax, val)
                    if set_dmin:
                        column.dmin = min(column.dmin, val)
            except KeyError:
                pass  # missing data - skip

        # Limit auto-generated scales with floor, ceiling and minrange.
        if column.ceiling is not None and column.max is None:
            column.dmax = min(column.dmax, float(column.ceiling))
        if column.floor is not None and column.min is None:
            column.dmin = max(column.dmin, float(column.floor))
        if column.minrange is not None:
            ddiff = column.dmax - column.dmin
            if ddiff < float(column.minrange):
                column.dmax = column.dmin + float(column.minrange)


def _collect_shared_keys(sections) -> Dict[str, Dict[str, Union[int, float]]]:
    # Collect settings for shared keys
    shared_keys: Dict[str, Dict[str, Union[int, float]]] = defaultdict(lambda: dict())
    for sec_idx, section in enumerate(sections):
        for _, column in section.column_by_key.items():
            sk: str = column.shared_key
            if sk is not None:
                shared_keys[sk]["dmax"] = max(
                    column.dmax,
                    shared_keys[sk].get("dmax", column.dmax),
                )
                shared_keys[sk]["dmin"] = min(
                    column.dmin,
                    shared_keys[sk].get("dmin", column.dmin),
                )
    return shared_keys
