"""
MultiQC datatable class, used by tables and violin plots
"""

import logging
import math
import re
from collections import defaultdict
from typing import Any, Callable, Dict, List, Mapping, NewType, Optional, Sequence, Set, Tuple, TypedDict, Union, cast

from pydantic import BaseModel, Field

from multiqc import config, report
from multiqc.plots.plotly.plot import PConfig
from multiqc.types import Anchor, ColumnKey, SampleGroup, SampleName
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


ColumnAnchor = NewType("ColumnAnchor", str)  # Unique within a table

ValueT = Union[int, float, str, bool]


def is_valid_value(val: Union[ValueT, None]) -> bool:
    """
    Run time check if the value is valid for a table cell. Duplicates the type hint, but
    used in case if a module ignores type checking.
    """
    return isinstance(val, (int, float, str, bool)) or val is None  # type: ignore


class ColumnDict(TypedDict, total=False):
    rid: ColumnAnchor  # namespace + short_rid = ID unique within a table
    title: str
    description: str
    scale: Union[str, bool]
    hidden: bool
    placement: float
    namespace: str
    color: Optional[str]
    colour: Optional[str]  # deprecated
    max: Optional[float]
    dmax: Optional[float]
    min: Optional[float]
    dmin: Optional[float]
    ceiling: Optional[float]
    floor: Optional[float]
    minrange: Optional[float]
    minRange: Optional[float]  # deprecated
    shared_key: Optional[str]
    tt_decimals: Optional[int]
    suffix: Optional[str]
    cond_formatting_colours: List[Dict[str, str]]
    cond_formatting_rules: Dict[str, List[Dict[str, Union[str, int, float]]]]
    bgcols: Dict[str, str]
    bars_zero_centrepoint: bool
    modify: Optional[Callable[[ValueT], ValueT]]
    format: Optional[Union[str, Callable[[ValueT], str]]]


class ColumnMeta(ValidatedConfig):
    """
    Column model class. Holds configuration for a single column in a table.
    """

    rid: ColumnAnchor  # namespace + short_rid = ID unique within a table
    title: str
    description: str
    scale: Union[str, bool]
    hidden: bool = False
    placement: float = 1000
    namespace: str = ""
    colour: Optional[str] = Field(None, deprecated="color")
    color: Optional[str] = None
    max: Optional[float] = None
    dmax: Optional[float] = None
    min: Optional[float] = None
    dmin: Optional[float] = None
    ceiling: Optional[float] = None
    floor: Optional[float] = None
    minRange: Optional[float] = Field(None, deprecated="minrange")
    minrange: Optional[float] = None
    shared_key: Optional[str] = None
    tt_decimals: Optional[int] = None
    suffix: Optional[str] = None
    cond_formatting_colours: List[Dict[str, str]] = []
    cond_formatting_rules: Dict[str, List[Dict[str, Union[str, int, float]]]] = {}
    bgcols: Dict[str, str] = {}
    bars_zero_centrepoint: bool = False
    modify: Optional[Callable[[ValueT], ValueT]] = None
    format: Optional[Union[str, Callable[[ValueT], str]]] = None

    @staticmethod
    def create(
        col_dict: ColumnDict,
        col_key: ColumnKey,  # to initialize rid
        sec_idx: int,  # to initialize the colour
        pconfig: TableConfig,  # plot config dictionary
        table_anchor: Anchor,
    ) -> "ColumnMeta":
        # Overwrite any header config if set in config
        if header_config := config.custom_table_header_config.get(pconfig.id, {}):
            if col_config := header_config.get(col_key, {}):
                for custom_k, custom_v in col_config.items():
                    col_dict[custom_k] = custom_v  # type: ignore

        namespace = col_dict.get("namespace", pconfig.namespace) or ""
        assert isinstance(namespace, str)

        unclean_rid = col_dict.get("rid") or col_key
        legacy_short_rid = re.sub(r"\W+", "_", str(unclean_rid)).strip().strip("_")  # User configs can still use it
        _rid = legacy_short_rid
        # Prefixing with namepsace to get a unique column ID within a table across all sections
        if namespace:
            ns_slugified = re.sub(r"\W+", "_", str(namespace)).strip().strip("_").lower()
            _rid = f"{ns_slugified}-{_rid}"
        col_dict["rid"] = ColumnAnchor(report.save_htmlid(_rid, scope=table_anchor))

        # Additionally override the header config assuming rid is used (for legacy docs)
        if header_config:
            if col_config := header_config.get(col_dict["rid"], {}):
                for custom_k, custom_v in col_config.items():
                    col_dict[custom_k] = custom_v  # type: ignore

        # Applying defaults presets for data keys if shared_key is set to base_count or read_count
        shared_key = col_dict.get("shared_key", None)
        if shared_key in ["read_count", "long_read_count", "base_count"]:
            shared_key_suffix = None
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
            if col_dict.get("modify") is None:
                col_dict["modify"] = lambda x: x * multiplier if isinstance(x, (int, float)) else x  # type: ignore  # noqa: E731
            if col_dict.get("min") is None:
                col_dict["min"] = 0
            if col_dict.get("format") is None and multiplier == 1:
                col_dict["format"] = "{:,d}"
            if col_dict.get("suffix") is None and shared_key_suffix is not None:
                col_dict["suffix"] = " " + shared_key_suffix

        col_dict.setdefault("min", pconfig.min)
        col_dict.setdefault("scale", pconfig.scale)
        col_dict.setdefault("description", col_dict.setdefault("title", str(col_key)))
        col_dict.setdefault("placement", 1000)

        # Overwrite (2nd time) any given config with table-level user config
        # This is to override column-specific values set by modules
        if pconfig.id in config.custom_plot_config:
            for cpc_k, cpc_v in config.custom_plot_config[pconfig.id].items():
                if isinstance(cpc_k, str) and cpc_k in ColumnMeta.model_fields.keys():
                    col_dict[cpc_k] = cpc_v  # type: ignore

        col: ColumnMeta = ColumnMeta(**col_dict)

        if col.color is None:
            cidx = sec_idx
            while cidx >= len(SECTION_COLORS):
                cidx -= len(SECTION_COLORS)
            col.color = SECTION_COLORS[cidx]

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
                        if _col_match(item_id2) and isinstance(visible, bool):
                            # Config has True = visible, False = Hidden. Here we're setting "hidden" which is inverse
                            col.hidden = not visible

            # Case-insensitive check if the outer key is a column ID
            elif _col_match(item_id) and isinstance(visibility, bool):
                # Config has True = visible, False = Hidden. Here we're setting "hidden" which is inverse
                col.hidden = not visibility

        # Also overwrite placement if set in config
        for item_id, item in config.table_columns_placement.items():
            if _ns_match(item_id) and isinstance(item, dict):
                for item_id2, placement in item.items():
                    if _col_match(item_id2) and isinstance(placement, (float, int)):
                        col.placement = float(placement)
            elif _col_match(item_id) and isinstance(item, (float, int)):
                col.placement = float(item)

        return col


class InputRow(BaseModel):
    """
    Row class. Holds configuration for a single row in a table (can be multiple for one sample)
    """

    sample: SampleName
    data: Dict[ColumnKey, Optional[ValueT]] = dict()

    def __init__(self, sample: SampleName, data: Mapping[Union[str, ColumnKey], Any]):
        super().__init__(
            sample=sample,
            data={ColumnKey(k): v for k, v in data.items() if is_valid_value(v)},
        )


ColumnKeyT = Union[str, ColumnKey]
GroupKeyT = Union[str, SampleGroup]
GroupT = Union[Mapping[ColumnKeyT, Optional[ValueT]], InputRow, Sequence[InputRow]]
SectionT = Mapping[GroupKeyT, GroupT]


class Row(BaseModel):
    """
    Processed row class. Holds configuration for a single row in a table (can be multiple for one sample).
    Contains raw, optionally modified, non-null values, and corresponding formatted values to display.
    """

    sample: SampleName
    # rows with original, unformatted values coming from modules:
    raw_data: Dict[ColumnKey, ValueT] = dict()
    # formatted rows (i.e. values are HTML strings to display in the table):
    formatted_data: Dict[ColumnKey, str] = dict()


class TableSection(BaseModel):
    """
    Table section class. Holds configuration for a single section in a table.
    """

    column_by_key: Dict[ColumnKey, ColumnMeta]
    rows_by_sgroup: Dict[SampleGroup, List[Row]] = defaultdict(list)


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

col_anchors_by_table: Dict[Anchor, Set[ColumnAnchor]] = defaultdict(set)


class DataTable(BaseModel):
    """
    Data table class. Prepares and holds data and configuration
    for either a table or a violin plot.
    """

    id: str
    anchor: Anchor
    pconfig: TableConfig

    sections: List[TableSection]
    headers_in_order: Dict[float, List[Tuple[int, ColumnKey]]]

    @staticmethod
    def create(
        data: Union[SectionT, List[SectionT]],
        table_id: str,
        table_anchor: Anchor,
        pconfig: TableConfig,
        headers: Optional[
            Union[
                List[Dict[ColumnKey, ColumnDict]],
                List[Dict[str, ColumnDict]],
                Dict[ColumnKey, ColumnDict],
                Dict[str, ColumnDict],
            ]
        ] = None,
    ) -> "DataTable":
        """Prepare data for use in a table or plot"""
        if headers is None:
            headers = []

        # Given one dataset - turn it into a list
        list_of_headers = cast(
            Union[List[Dict[ColumnKey, ColumnDict]], List[Dict[str, ColumnDict]]],
            headers if isinstance(headers, list) else [headers],
        )
        del headers

        # Each section to have a list of groups (even if there is just one element in a group)
        input_section: SectionT
        input_group: GroupT
        unified_sections__with_nulls: List[Dict[SampleGroup, List[InputRow]]] = []
        for input_section in data if isinstance(data, list) else [data]:
            rows_by_group: Dict[SampleGroup, List[InputRow]] = {}
            for g_name, input_group in input_section.items():
                g_name = SampleGroup(str(g_name))  # Make sure sample names are strings
                if isinstance(input_group, dict):  # just one row, defined as a mapping from metric to value
                    # Remove non-scalar values for table cells
                    rows_by_group[g_name] = [InputRow(sample=SampleName(g_name), data=input_group)]
                elif isinstance(input_group, list):  # multiple rows, each defined as a mapping from metric to value
                    rows_by_group[g_name] = input_group
                else:
                    assert isinstance(input_group, InputRow)
                    rows_by_group[g_name] = [input_group]
            unified_sections__with_nulls.append(rows_by_group)
        del data

        # Go through each table section and create a list of Section objects
        sections: List[TableSection] = []
        for sec_idx, rows_by_sname__with_nulls in enumerate(unified_sections__with_nulls):
            header_by_key: Union[Dict[ColumnKey, ColumnDict], Dict[str, ColumnDict]] = (
                list_of_headers[sec_idx] if sec_idx < len(list_of_headers) else dict()
            )
            if not header_by_key:
                pconfig.only_defined_headers = False

            column_by_key: Dict[ColumnKey, ColumnMeta] = dict()
            col_dict_by_key_copy: Dict[ColumnKey, ColumnDict] = _get_or_create_headers(
                rows_by_sname__with_nulls, header_by_key, pconfig
            )
            for col_key, col_dict in col_dict_by_key_copy.items():
                column_by_key[col_key] = ColumnMeta.create(
                    col_dict=col_dict, col_key=col_key, sec_idx=sec_idx, pconfig=pconfig, table_anchor=table_anchor
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
        headers_in_order: Dict[float, List[Tuple[int, ColumnKey]]] = defaultdict(list)
        for sec_idx, section in enumerate(sections):
            for col_key, column in section.column_by_key.items():
                if column.shared_key is not None:
                    column.dmax = shared_keys[column.shared_key].get("dmax")
                    column.dmin = shared_keys[column.shared_key].get("dmin")

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

    def get_headers_in_order(self) -> List[Tuple[int, ColumnKey, ColumnMeta]]:
        """
        Gets the headers in the order they want to be displayed.
        Returns a list of triplets: (bucket_idx, key, header_info)
        """
        res: List[Tuple[int, ColumnKey, ColumnMeta]] = list()
        # Scan through self.headers_in_order and just bolt on the actual header info
        placement: float
        for placement in sorted(self.headers_in_order.keys()):
            for section_idx, col_key in self.headers_in_order[placement]:
                res.append((section_idx, col_key, self.sections[section_idx].column_by_key[col_key]))
        return res


def _get_or_create_headers(
    rows_by_sample: Dict[SampleGroup, List[InputRow]],
    header_by_key: Union[Mapping[str, ColumnDict], Mapping[ColumnKey, ColumnDict]],
    pconfig: TableConfig,
) -> Dict[ColumnKey, ColumnDict]:
    """
    Process and populate headers, if missing or incomplete.
    """
    # Make a copy to keep the input immutable.
    header_by_key_copy = {ColumnKey(k): h for k, h in header_by_key.items()}
    if not pconfig.only_defined_headers:
        # Get additional header keys from the data
        col_ids: List[ColumnKey] = list(header_by_key_copy.keys())
        # Get the keys from the data
        for _, rows in rows_by_sample.items():
            for row in rows:
                for col_id in row.data.keys():
                    if col_id not in col_ids:
                        col_ids.append(col_id)

        # Create empty header configs for each new data key
        for col_id in col_ids:
            if col_id not in header_by_key:
                header_by_key_copy[col_id] = {}

    # Check that we have some data in each column
    empties: List[ColumnKey] = list()
    for col_id in header_by_key_copy.keys():
        n = 0
        for _, rows in rows_by_sample.items():
            for row in rows:
                if col_id in row.data.keys():
                    n += 1

        if n == 0:
            empties.append(col_id)

    # Remove empty columns
    for empty_col_id in empties:
        logger.debug(
            f"Table key '{empty_col_id}' not found in data for '{pconfig.id}'. Skipping. Check for possible typos "
            f"between data keys and header keys"
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
            logger.error(f"Error modifying table value '{column.rid}': '{val}'. {e}")
    # section.raw_data[s_name][col_key] = val

    # Now also calculate formatted values
    valstr = str(val)
    fmt: Union[None, str, Callable[[ValueT], str]] = column.format
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
    col_key: ColumnKey,
    rows_by_sample: Dict[SampleGroup, List[Row]],
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
        for _, rows in rows_by_sample.items():
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


def _collect_shared_keys(sections: List[TableSection]) -> Dict[str, Dict[str, Union[int, float]]]:
    # Collect settings for shared keys
    shared_keys: Dict[str, Dict[str, Union[int, float]]] = defaultdict(lambda: dict())
    for _, section in enumerate(sections):
        for _, column in section.column_by_key.items():
            sk: Optional[str] = column.shared_key
            if sk is not None:
                sk_dmax: Optional[float] = shared_keys[sk].get("dmax")
                if sk_dmax is not None and column.dmax is not None:
                    shared_keys[sk]["dmax"] = max(column.dmax, sk_dmax)
                elif sk_dmax is None and column.dmax is not None:
                    shared_keys[sk]["dmax"] = column.dmax
                else:
                    pass

                sk_dmin: Optional[float] = shared_keys[sk].get("dmin")
                if sk_dmin is not None and column.dmin is not None:
                    shared_keys[sk]["dmin"] = min(column.dmin, sk_dmin)
                elif sk_dmin is None and column.dmin is not None:
                    shared_keys[sk]["dmin"] = column.dmin
                else:
                    pass

    return shared_keys
