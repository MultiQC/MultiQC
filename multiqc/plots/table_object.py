"""
MultiQC datatable class, used by tables and violin plots
"""

import logging
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Mapping, NewType, Optional, Sequence, Set, Tuple, TypedDict, Union, cast
from pathlib import Path

from natsort import natsorted
from pydantic import BaseModel, Field

from multiqc import config, report
from multiqc.plots.plot import PConfig
from multiqc.types import Anchor, ColumnKey, SampleGroup, SampleName, SectionKey
from multiqc.utils import mqc_colour
from multiqc.utils.material_icons import get_material_icon
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
    parse_numeric: bool = True
    rows_are_samples: bool = True
    flat_if_very_large: bool = False

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        super().__init__(path_in_cfg=path_in_cfg or ("table",), **data)


ColumnAnchor = NewType("ColumnAnchor", str)  # Unique within a table


ValueT = Union[int, float, str, bool]


@dataclass
class Cell:
    raw: ValueT
    mod: ValueT
    fmt: str


ExtValueT = Union[int, float, str, bool, Cell]


def is_valid_value(val: Union[ExtValueT, None]) -> bool:
    """
    Run time check if the value is valid for a table cell. Duplicates the type hint, but
    used in case if a module ignores type checking.
    """
    return isinstance(val, (int, float, str, bool, Cell)) or val is None  # type: ignore


class ColumnDict(TypedDict, total=False):
    rid: ColumnAnchor  # namespace + short_rid = ID unique within a table
    clean_rid: ColumnAnchor  # can differ when rid is provided by user
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
    clean_rid: ColumnAnchor  # can differ when rid is provided by user
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
        header_config = None
        for _id in [table_anchor, pconfig.id]:  # Back-compatibility using pconfig.id instead of table_anchor
            if _header_config := config.custom_table_header_config.get(_id, {}):
                header_config = _header_config
                if col_config := header_config.get(col_key, {}):
                    for custom_k, custom_v in col_config.items():
                        col_dict[custom_k] = custom_v  # type: ignore
                break

        namespace = col_dict.get("namespace", pconfig.namespace) or ""
        assert isinstance(namespace, str)

        unclean_rid = col_dict.get("rid") or col_key
        legacy_short_rid = re.sub(r"\W+", "_", str(unclean_rid)).strip().strip("_")  # User configs can still use it
        if "clean_rid" not in col_dict:
            _rid = legacy_short_rid
            # Prefixing with namepsace to get a unique column ID within a table across all sections
            if namespace:
                ns_slugified = re.sub(r"\W+", "_", str(namespace)).strip().strip("_").lower()
                _rid = f"{ns_slugified}-{_rid}"
            col_dict["clean_rid"] = ColumnAnchor(report.save_htmlid(_rid, scope=table_anchor))
        if "rid" not in col_dict:
            col_dict["rid"] = col_dict["clean_rid"]

        # Additionally override the header config assuming rid is used (for legacy docs)
        if header_config:
            if col_config := header_config.get(col_dict["clean_rid"], {}):
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
        for _id in [table_anchor, pconfig.id]:  # Back-compatibility using pconfig.id instead of table_anchor
            if config_dict := config.custom_plot_config.get(_id):
                for cpc_k, cpc_v in config_dict.items():
                    if isinstance(cpc_k, str) and cpc_k in ColumnMeta.model_fields.keys():
                        col_dict[cpc_k] = cpc_v  # type: ignore
                break

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
                    table_anchor,
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
                    col.clean_rid,
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

    def __init__(self, path_in_cfg: Optional[Tuple[str, ...]] = None, **data):
        super().__init__(path_in_cfg=path_in_cfg or ("table", "column"), **data)


class InputRow(BaseModel):
    """
    Row class. Holds configuration for a single row in a table (can be multiple for one sample)
    """

    sample: SampleName
    data: Dict[ColumnKey, Optional[ExtValueT]] = Field(default_factory=dict)

    def __init__(self, sample: SampleName, data: Mapping[Union[str, ColumnKey], Any]):
        super().__init__(
            sample=sample,
            data={ColumnKey(k): v for k, v in data.items() if is_valid_value(v)},
        )


ColumnKeyT = Union[str, ColumnKey]
GroupKeyT = Union[str, SampleGroup]
GroupT = Union[Mapping[ColumnKeyT, Optional[ExtValueT]], InputRow, Sequence[InputRow]]
SectionT = Mapping[GroupKeyT, GroupT]


class Row(BaseModel):
    """
    Processed row class. Holds configuration for a single row in a table (can be multiple for one sample).
    Contains raw, optionally modified, non-null values, and corresponding formatted values to display.
    """

    sample: SampleName
    data: Dict[ColumnKey, Cell] = dict()


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

    section_by_id: Dict[SectionKey, TableSection]
    headers_in_order: Dict[float, List[Tuple[int, ColumnKey]]]

    def is_empty(self) -> bool:
        """
        Check if the table is empty.
        """
        return not self.section_by_id or not any(
            section.rows_by_sgroup.values() for section in self.section_by_id.values()
        )

    @staticmethod
    def create(
        data: Dict[SectionKey, SectionT],
        table_id: str,
        table_anchor: Anchor,
        pconfig: TableConfig,
        headers: Dict[SectionKey, Dict[ColumnKey, ColumnDict]],
    ) -> "DataTable":
        """Prepare data for use in a table or plot"""
        # Violin plot's PConfig does this for the plot ID. We also do that second turn time for
        # the table anchor because that's the ID that is shown in the Configure Columns modal
        if table_anchor in config.custom_plot_config:
            for k, v in config.custom_plot_config[table_anchor].items():
                if isinstance(k, str) and k in pconfig.model_fields:
                    setattr(pconfig, k, v)

        # Each section to have a list of groups (even if there is just one element in a group)
        input_section_key: SectionKey
        input_section: SectionT
        input_group: GroupT
        unified_sections__with_nulls: Dict[SectionKey, Dict[SampleGroup, List[InputRow]]] = {}
        for input_section_key, input_section in data.items():
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
            unified_sections__with_nulls[input_section_key] = rows_by_group
        del data

        # Go through each table section and create a list of Section objects
        sections: Dict[SectionKey, TableSection] = {}
        for sec_idx, (section_key, rows_by_sname__with_nulls) in enumerate(unified_sections__with_nulls.items()):
            header_by_key: Union[Dict[ColumnKey, ColumnDict], Dict[str, ColumnDict]] = (
                headers.get(SectionKey(section_key)) or dict()
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
                        val_obj = _process_and_format_value(
                            optional_val, column_by_key[col_key], parse_numeric=pconfig.parse_numeric
                        )
                        row.data[col_key] = val_obj
                    if row.data:
                        section.rows_by_sgroup[g_name].append(row)

            # Remove empty groups:
            section.rows_by_sgroup = {sname: rows for sname, rows in section.rows_by_sgroup.items() if rows}

            # Work out max and min value if not given:
            for col_key, column in column_by_key.items():
                _determine_dmin_and_dmax(column, col_key, section.rows_by_sgroup)

            sections[section_key] = section

        del unified_sections__with_nulls

        shared_keys: Dict[str, Dict[str, Union[int, float]]] = _collect_shared_keys(sections)

        # Overwrite shared key settings and at the same time assign to buckets for sorting
        # So the final ordering is:
        #   placement > section > explicit_ordering
        # Of course, the user can shuffle these manually.
        headers_in_order: Dict[float, List[Tuple[int, ColumnKey]]] = defaultdict(list)
        for sec_idx, section in enumerate(sections.values()):
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
        for section in sections.values():
            for column in section.column_by_key.values():
                column.modify = None
                if not isinstance(column.format, str):
                    column.format = None

        # Assign to class
        return DataTable(
            id=table_id,
            anchor=table_anchor,
            section_by_id=sections,
            headers_in_order=headers_in_order,
            pconfig=pconfig,
        )

    def get_headers_in_order(self, keep_hidden: bool = True) -> List[Tuple[int, ColumnKey, ColumnMeta]]:
        """
        Gets the headers in the order they want to be displayed.
        Returns a list of triplets: (bucket_idx, key, header_info)
        """
        res: List[Tuple[int, ColumnKey, ColumnMeta]] = list()
        # Scan through self.headers_in_order and just bolt on the actual header info
        placement: float
        for placement in sorted(self.headers_in_order.keys()):
            for section_idx, col_key in self.headers_in_order[placement]:
                if keep_hidden or not list(self.section_by_id.values())[section_idx].column_by_key[col_key].hidden:
                    res.append(
                        (section_idx, col_key, list(self.section_by_id.values())[section_idx].column_by_key[col_key])
                    )
        return res

    def merge(self, new_dt: "DataTable"):
        """
        Extend the existing DataTable with new data.

        Sections are matched by their column keys rather than their position in the list.
        This allows for merging sections even if they appear in different order.

        New samples and columns are added at the end, while existing sample/column combinations
        are overwritten with new data.
        """

        # Find matching sections and merge them
        for section_id, new_section in new_dt.section_by_id.items():
            if section_id in self.section_by_id:
                # Merge matching sections
                section = self.section_by_id[section_id]

                # Update column metadata
                for col_key, new_col_meta in new_section.column_by_key.items():
                    if col_key in section.column_by_key:
                        # Update existing column metadata fields
                        existing_col = section.column_by_key[col_key]
                        for field_name, field_value in new_col_meta.model_dump().items():
                            setattr(existing_col, field_name, field_value)
                    else:
                        # Add new column
                        section.column_by_key[col_key] = new_col_meta

                # Update row data
                for group_name, new_rows in new_section.rows_by_sgroup.items():
                    existing_rows = section.rows_by_sgroup.get(group_name, [])
                    existing_samples = {row.sample for row in existing_rows}

                    for new_row in new_rows:
                        if new_row.sample in existing_samples:
                            # Update existing sample data
                            for i, row in enumerate(existing_rows):
                                if row.sample == new_row.sample:
                                    row.data.update(new_row.data)
                                    existing_rows[i] = row
                                    break
                        else:
                            # Add new sample
                            existing_rows.append(new_row)

                    section.rows_by_sgroup[group_name] = existing_rows

                # Work out max and min value if not given:
                for col_key, column in section.column_by_key.items():
                    _determine_dmin_and_dmax(column, col_key, section.rows_by_sgroup)

            else:
                # Add new section if no match found
                self.section_by_id[section_id] = new_section

        # Rebuild headers_in_order based on updated columns
        self._rebuild_headers_in_order()

    def _rebuild_headers_in_order(self) -> None:
        """
        Rebuild the headers_in_order dictionary based on current sections and columns.
        """
        self.headers_in_order = defaultdict(list)
        for sec_idx, section in enumerate(self.section_by_id.values()):
            for col_key, column in section.column_by_key.items():
                self.headers_in_order[column.placement].append((sec_idx, col_key))


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


def _process_and_format_value(val: ExtValueT, column: ColumnMeta, parse_numeric: bool = True) -> Cell:
    """
    Takes row value, applies "modify" functions and "format" string, and returns a tuple:
    the modified value and its formatted string.

    "parse_numeric=False" assumes that the numeric values are already pre-parsed
    """
    if isinstance(val, Cell):  # aready formatted
        return val

    # Try parse as a number
    if parse_numeric:
        if str(val).isdigit():
            val = int(val)
        else:
            try:
                val = float(val)
            except ValueError:
                pass

    val_unmodified = val
    # Apply modify
    if column.modify:
        # noinspection PyBroadException
        try:
            val = column.modify(val)
        except Exception as e:  # User-provided modify function can raise any exception
            logger.error(f"Error modifying table value '{column.clean_rid}': '{val}'. {e}")

    # Values can be quoted to avoid parsing as a number, so need to remove those quotes
    if isinstance(val, str) and (
        val.startswith('"') and val.endswith('"') or val.startswith("'") and val.endswith("'")
    ):
        val = val[1:-1]

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
                logger.debug(f"Error applying format to table value '{column.clean_rid}': '{val}'. {e}")
        elif isinstance(val, (int, float)):
            if fmt == r"{:,d}":
                val = round(val)
            try:
                # If format is decimal and value is float, try rounding to int first
                if isinstance(val, float) and "d" in fmt:
                    try:
                        val = int(round(val))
                    except (ValueError, OverflowError):
                        pass
                valstr = fmt.format(val)
            except Exception as e:
                logger.debug(
                    f"Error applying format string '{fmt}' to table value '{column.clean_rid}': '{val}'. {e}. "
                    f"Check if your format string is correct."
                )
    return Cell(raw=val_unmodified, mod=val, fmt=valstr)


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
            v_by_col = rows[0].data
            try:
                val = v_by_col[col_key].mod if v_by_col[col_key].mod is not None else v_by_col[col_key].raw
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


def _collect_shared_keys(sections: Dict[SectionKey, TableSection]) -> Dict[str, Dict[str, Union[int, float]]]:
    # Collect settings for shared keys
    shared_keys: Dict[str, Dict[str, Union[int, float]]] = defaultdict(lambda: dict())
    for _, section in sections.items():
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


def render_html(
    dt: DataTable,
    violin_anchor: Anchor,
    module_anchor: Anchor,
    section_anchor: Anchor,
    add_control_panel: bool = True,
) -> Tuple[str, str]:
    """
    Build HTML for a MultiQC table, and HTML for the modal for configuring the table.
    :param dt: MultiQC datatable object
    :param violin_anchor: optional, will add a button to switch to a violin plot with this ID
    :param add_control_panel: whether to add the control panel with buttons above the table
    """

    col_to_th: Dict[ColumnAnchor, str] = dict()
    col_to_modal_headers: Dict[ColumnAnchor, str] = dict()
    col_to_hidden: Dict[ColumnAnchor, bool] = dict()
    group_to_sample_to_anchor_to_td: Dict[SampleGroup, Dict[SampleName, Dict[ColumnAnchor, str]]] = defaultdict(
        lambda: defaultdict(dict)
    )
    group_to_sample_to_anchor_to_val: Dict[SampleGroup, Dict[SampleName, Dict[ColumnAnchor, ValueT]]] = defaultdict(
        lambda: defaultdict(dict)
    )
    group_to_sample_to_nice_name_to_val: Dict[SampleGroup, Dict[SampleName, Dict[str, ValueT]]] = defaultdict(
        lambda: defaultdict(dict)
    )
    group_to_sorting_to_anchor_to_val: Dict[SampleGroup, Dict[ColumnAnchor, ValueT]] = defaultdict(dict)
    group_to_sample_to_anchor_to_empty: Dict[SampleGroup, Dict[SampleName, Dict[ColumnAnchor, bool]]] = defaultdict(
        lambda: defaultdict(dict)
    )
    # empty_cells: Dict[ColumnKeyT, str] = dict()
    hidden_cols = 1
    table_title = dt.pconfig.title

    def escape(s: str) -> str:
        return s.replace('"', "&quot;").replace("'", "&#39;").replace("<", "&lt;").replace(">", "&gt;")

    for idx, col_key, header in dt.get_headers_in_order():
        col_anchor: ColumnAnchor = header.clean_rid

        # Build the table header cell
        shared_key = ""
        if header.shared_key is not None:
            shared_key = f" data-shared-key={header.shared_key}"

        td_hide_cls = ""
        tr_muted_cls = ""
        checked = ' checked="checked"'
        if header.hidden:
            td_hide_cls = " column-hidden"
            tr_muted_cls = " text-muted"
            checked = ""
            hidden_cols += 1

        data_attr = (
            f'data-dmax="{header.dmax}" data-dmin="{header.dmin}" data-namespace="{header.namespace}" {shared_key}'
        )

        ns = f"{header.namespace}: " if header.namespace else ""
        cell_contents = (
            f'<span class="mqc_table_tooltip" title="{ns}{header.description}" data-html="true">{header.title}</span>'
        )

        col_to_th[col_anchor] = (
            f'<th id="header_{col_anchor}" class="{col_anchor}{td_hide_cls}" {data_attr}>{cell_contents}</th>'
        )
        col_to_hidden[col_anchor] = header.hidden

        # Build the modal table row
        data = f"data-table-anchor='{dt.anchor}'"
        if violin_anchor:
            data += f" data-violin-anchor='{violin_anchor}'"
        col_to_modal_headers[col_anchor] = f"""
        <tr class="{col_anchor}{tr_muted_cls}" style="background-color: rgba({header.color}, 0.15);">
          <td class="sorthandle ui-sortable-handle">||</span></td>
          <td style="text-align:center;">
            <input class="mqc_table_col_visible" type="checkbox" {checked} value="{col_anchor}" {data}>
          </td>
          <td>{header.namespace}</td>
          <td>{header.title}</td>
          <td>{header.description}</td>
          <td><code>{col_anchor}</code></td>
          <td>{header.shared_key or ""}</td>
        </tr>"""

        # Make a colour scale
        c_scale = None
        if isinstance(header.scale, str):
            c_scale = mqc_colour.mqc_colour_scale(
                name=header.scale,
                minval=header.dmin,
                maxval=header.dmax,
                id=dt.id,
            )

        # Collect conditional formatting config
        cond_formatting_rules: Dict[str, Dict[str, List[Dict[str, Union[str, int, float]]]]] = {}
        if header.cond_formatting_rules:
            cond_formatting_rules[col_anchor] = header.cond_formatting_rules
        cond_formatting_rules.update(config.table_cond_formatting_rules)

        cond_formatting_colours = header.cond_formatting_colours
        cond_formatting_colours.extend(config.table_cond_formatting_colours)

        # Add the data table cells
        section = list(dt.section_by_id.values())[idx]
        for group_name, group_rows in section.rows_by_sgroup.items():
            for row_idx, row in enumerate(group_rows):
                if col_key not in row.data:
                    continue

                val: ValueT = row.data[col_key].mod
                valstr: str = row.data[col_key].fmt

                group_to_sample_to_anchor_to_val[group_name][row.sample][col_anchor] = val
                group_to_sample_to_nice_name_to_val[group_name][row.sample][col_key] = val

                if c_scale and c_scale.name not in c_scale.qualitative_scales:
                    dmin = header.dmin
                    dmax = header.dmax
                    if dmin is not None and dmax is not None and dmax != dmin:
                        try:
                            val_float = float(val)
                        except ValueError:
                            percentage = 0.0
                        else:
                            percentage = ((val_float - dmin) / (dmax - dmin)) * 100
                            # Treat 0 as 0-width and make bars width of absolute value
                            if header.bars_zero_centrepoint:
                                dmax = max(abs(dmin), abs(dmax))
                                dmin = 0
                                percentage = ((abs(val_float) - dmin) / (dmax - dmin)) * 100
                            percentage = min(percentage, 100)
                            percentage = max(percentage, 0)
                    else:
                        percentage = 0.0
                else:
                    percentage = 100.0

                # This is horrible, but Python locale settings are worse
                if config.thousandsSep_format is None:
                    config.thousandsSep_format = '<span class="mqc_small_space"></span>'
                if config.decimalPoint_format is None:
                    config.decimalPoint_format = "."
                valstr = valstr.replace(".", "DECIMAL").replace(",", "THOUSAND")
                valstr = valstr.replace("DECIMAL", config.decimalPoint_format).replace(
                    "THOUSAND", config.thousandsSep_format
                )

                suffix = header.suffix
                if suffix:
                    # Add a space before the suffix, but not as an actual character, so ClipboardJS would copy
                    # the whole value without the space. Also, remove &nbsp; that we don't want ClipboardJS to copy.
                    suffix = suffix.replace("&nbsp;", " ").strip()
                    valstr += "<span class='mqc_small_space'></span>" + suffix

                # Conditional formatting
                # Build empty dict for color formatting matches:
                cmatches = {}
                for cfc in cond_formatting_colours:
                    for cfc_key in cfc:
                        cmatches[cfc_key] = False
                # Find general rules followed by column-specific rules
                for cfk in ["all_columns", str(col_anchor), str(dt.id)]:
                    if cfk in cond_formatting_rules:
                        # Loop through match types
                        for ftype in cmatches.keys():
                            # Loop through array of comparison types
                            for cmp in cond_formatting_rules[cfk].get(ftype, []):
                                try:
                                    # Each comparison should be a dict with single key: val
                                    if "s_eq" in cmp and str(cmp["s_eq"]).lower() == str(val).lower():
                                        cmatches[ftype] = True
                                    if "s_contains" in cmp and str(cmp["s_contains"]).lower() in str(val).lower():
                                        cmatches[ftype] = True
                                    if "s_ne" in cmp and str(cmp["s_ne"]).lower() != str(val).lower():
                                        cmatches[ftype] = True
                                    if "eq" in cmp and float(val) == float(cmp["eq"]):
                                        cmatches[ftype] = True
                                    if "ne" in cmp and float(val) != float(cmp["ne"]):
                                        cmatches[ftype] = True
                                    if "gt" in cmp and float(val) > float(cmp["gt"]):
                                        cmatches[ftype] = True
                                    if "lt" in cmp and float(val) < float(cmp["lt"]):
                                        cmatches[ftype] = True
                                    if "ge" in cmp and float(val) >= float(cmp["ge"]):
                                        cmatches[ftype] = True
                                    if "le" in cmp and float(val) <= float(cmp["le"]):
                                        cmatches[ftype] = True
                                except Exception as e:
                                    logger.warning(
                                        f"Not able to apply table conditional formatting to '{val}' ({cmp}): {e}"
                                    )
                # Apply HTML in order of config keys
                badge_col = None
                for cfc in cond_formatting_colours:
                    for cfc_key in cfc:  # should always be one, but you never know
                        if cmatches[cfc_key]:
                            badge_col = cfc[cfc_key]
                if badge_col is not None:
                    valstr = f'<span class="badge" style="background-color:{badge_col}">{valstr}</span>'

                # Determine background color based on scale. Only relevant for hashable values. If value is for some
                # reason a dict or a list, it's not hashable and the logic determining the color will not work.
                hashable = True
                try:
                    hash(val)
                except TypeError:
                    hashable = False
                    logger.warning(
                        f"Value {val} is not hashable for table {dt.anchor}, column {col_key}, sample {row.sample}"
                    )

                sorting_val = group_to_sorting_to_anchor_to_val.get(group_name, {}).get(col_anchor)
                if sorting_val is None:
                    group_to_sorting_to_anchor_to_val[group_name][col_anchor] = val
                    sorting_val = val

                # Categorical background colours supplied
                if isinstance(val, str) and val in header.bgcols.keys():
                    col = f'style="background-color:{header.bgcols[val]} !important;"'
                    group_to_sample_to_anchor_to_td[group_name][row.sample][col_anchor] = (
                        f'<td data-sorting-val="{escape(str(sorting_val))}" class="{col_anchor} {td_hide_cls}" {col}>{valstr}</td>'
                    )

                # Build table cell background colour bar
                elif hashable and header.scale:
                    if c_scale is not None:
                        col = " background-color:{} !important;".format(
                            c_scale.get_colour(val, source=f'Table "{dt.anchor}", column "{col_key}"')
                        )
                    else:
                        col = ""
                    bar_html = f'<span class="bar" style="width:{percentage}%;{col}"></span>'
                    val_html = f'<span class="val">{valstr}</span>'
                    wrapper_html = f'<div class="wrapper">{bar_html}{val_html}</div>'

                    group_to_sample_to_anchor_to_td[group_name][row.sample][col_anchor] = (
                        f'<td data-sorting-val="{escape(str(sorting_val))}" class="data-coloured {col_anchor} {td_hide_cls}">{wrapper_html}</td>'
                    )

                # Scale / background colours are disabled
                else:
                    group_to_sample_to_anchor_to_td[group_name][row.sample][col_anchor] = (
                        f'<td data-sorting-val="{escape(str(sorting_val))}" class="{col_anchor} {td_hide_cls}">{valstr}</td>'
                    )

                # Is this cell hidden or empty?
                group_to_sample_to_anchor_to_empty[group_name][row.sample][col_anchor] = (
                    header.hidden or str(val).strip() == ""
                )

        # Remove header if we don't have any filled cells for it
        sum_vals = 0
        for g, rows_by_sample in group_to_sample_to_anchor_to_td.items():
            sum_vals += sum([len(rows) for rows in rows_by_sample.values()])
        if sum_vals == 0:
            if header.hidden:
                hidden_cols -= 1
            col_to_th.pop(col_anchor, None)
            col_to_modal_headers.pop(col_anchor, None)
            logger.debug(f"Removing header {col_key} from table, as no data")

    # Put everything together
    html = ""

    # Buttons above the table
    if not config.simple_output and add_control_panel:
        # Copy Table Button
        buttons: List[str] = []

        buttons.append(
            f"""
        <button
            type="button"
            class="mqc_table_copy_btn btn btn-outline-secondary btn-sm"
            data-clipboard-target="table#{dt.anchor}"
            data-bs-toggle="tooltip"
            title="Copy table into clipboard suitable to be pasted into Excel or Google Sheets"
        >
            {get_material_icon("mdi:content-copy", 16)} Copy table
        </button>
        """
        )

        # Configure Columns Button
        if len(col_to_th) > 1:
            # performance degrades substantially when configuring thousands of columns
            # it is effectively unusable.
            disabled_class = ""
            disabled_attrs = ""
            if _is_configure_columns_disabled(len(col_to_th)):
                disabled_class = "mqc_table_tooltip"
                disabled_attrs = 'disabled title="Table is too large to configure columns"'

            buttons.append(
                f"""
            <button type="button" class="mqc_table_config_modal_btn btn btn-outline-secondary btn-sm {disabled_class}" data-bs-toggle="modal"
                data-bs-target="#{dt.anchor}_config_modal" {disabled_attrs} title="Configure visibility and ordering of columns">
                                    {get_material_icon("mdi:view-column", 16)} Configure columns
            </button>
            """
            )

        # Sort By Highlight button
        buttons.append(
            f"""
            <button type="button" class="mqc_table_sortHighlight btn btn-outline-secondary btn-sm"
                data-table-anchor="{dt.anchor}" data-direction="desc" style="display:none;" data-bs-toggle="tooltip" title="Place highlighted samples on top">
                                        {get_material_icon("mdi:sort", 16)} Sort by highlight
            </button>
        """
        )

        # Scatter Plot Button
        if len(col_to_th) > 1:
            buttons.append(
                f"""
                <button type="button" class="mqc_table_make_scatter btn btn-outline-secondary btn-sm"
                data-bs-toggle="modal" data-bs-target="#table_scatter_modal" data-table-anchor="{dt.anchor}" title="Visualize pairs of values on a scatter plot">
                                            {get_material_icon("mdi:chart-scatter-plot", 16)} Scatter plot
                </button>
                """
            )

        if violin_anchor is not None:
            buttons.append(
                f"""
                <button type="button" class="mqc-table-to-violin btn btn-outline-secondary btn-sm"
                data-table-anchor="{dt.anchor}" data-violin-anchor="{violin_anchor}" data-bs-toggle="tooltip" title="View as a violin plot">
                                            {get_material_icon("mdi:violin", 16)} Violin plot
                </button>
                """
            )

        buttons.append(
            f"""
        <button type="button" class="export-plot btn btn-outline-secondary btn-sm"
            data-plot-anchor="{violin_anchor or dt.anchor}" data-type="table" data-bs-toggle="tooltip" title="Show export options"
        >Export as CSV...</button>
        """
        )

        # "Showing x of y columns" text
        not_empty_rows_bool_vector = [
            all(group_to_sample_to_anchor_to_empty[s_name].values()) for s_name in group_to_sample_to_anchor_to_empty
        ]
        n_visible_rows = len([x for x in not_empty_rows_bool_vector if x is True])

        # Visible rows
        t_showing_rows_txt = f'Showing <sup id="{dt.anchor}_numrows" class="mqc_table_numrows">{n_visible_rows}</sup>/<sub>{len(group_to_sample_to_anchor_to_td)}</sub> rows'

        # How many columns are visible?
        ncols_vis = (len(col_to_th) + 1) - hidden_cols
        t_showing_cols_txt = ""
        if len(col_to_th) > 1:
            t_showing_cols_txt = f' and <sup id="{dt.anchor}_numcols" class="mqc_table_numcols">{ncols_vis}</sup>/<sub>{len(col_to_th)}</sub> columns'

        # Build table header text
        buttons.append(
            f"""
        <small id="{dt.anchor}_numrows_text" class="mqc_table_numrows_text">{t_showing_rows_txt}{t_showing_cols_txt}.</small>
        """
        )

        if not config.no_ai:
            seqera_ai_icon = (
                Path(__file__).parent.parent / "templates/default/assets/img/Seqera_AI_icon.svg"
            ).read_text()
            buttons.append(
                f"""
            <div class="ai-plot-buttons-container" style="float: right">
                <button
                    class="btn btn-outline-secondary btn-sm ai-copy-content ai-copy-content-table ai-copy-button-wrapper"
                    data-section-anchor="{section_anchor}"
                    data-plot-anchor="{violin_anchor}"
                    data-module-anchor="{module_anchor}"
                    data-plot-view="table"
                    type="button"
                    data-bs-toggle="tooltip"
                    title="Copy table data for use with AI tools like ChatGPT"
                >
                    {seqera_ai_icon}
                    <span class="button-text">Copy Prompt</span>
                </button>
                <button
                    class="btn btn-outline-secondary btn-sm ai-generate-button ai-generate-button-table ai-generate-button-wrapper"
                    data-response-div="{section_anchor}_ai_summary_response"
                    data-error-div="{section_anchor}_ai_summary_error"
                    data-disclaimer-div="{section_anchor}_ai_summary_disclaimer"
                    data-continue-in-chat-button="{section_anchor}_ai_summary_continue_in_chat"
                    data-detailed-analysis-div="{section_anchor}_ai_summary_detailed_analysis_response"
                    data-wrapper-div="{section_anchor}_ai_summary_wrapper"
                    data-section-anchor="{section_anchor}"
                    data-plot-anchor="{violin_anchor}"
                    data-module-anchor="{module_anchor}"
                    data-plot-view="table"
                    data-action="generate"
                    data-clear-text="Clear summary"
                    type="button"
                    data-bs-toggle="tooltip"
                    aria-controls="{dt.anchor}_ai_summary_wrapper"
                    title="Dynamically generate AI summary for this table"
                >
                    {seqera_ai_icon}
                    <span class="button-text">Summarize table</span>
                </button>
            </div>
            """
            )

        panel = "\n".join(buttons)
        html += f"""
        <div class='row mqc_table_control_buttons'>\n<div class='col-12'>\n{panel}\n</div>\n</div>
        """

    # Build the table itself
    collapse_class = (
        "mqc-table-collapse" if len(group_to_sample_to_anchor_to_td) > 10 and config.collapse_tables else ""
    )
    html += f"""
        <div id="{dt.anchor}_container" class="mqc_table_container">
            <div class="table-responsive mqc-table-responsive {collapse_class}" data-collapsed="{str(collapse_class != "").lower()}">
                <table id="{dt.anchor}" class="table table-sm mqc_table mqc_per_sample_table" data-title="{table_title}" data-sortlist="{_get_sortlist_js(dt)}">
        """

    # Build the header row
    col1_header = dt.pconfig.col1_header
    html += f'<thead><tr><th class="rowheader">{col1_header}</th>{"".join(col_to_th.values())}</tr></thead>'

    # Build the table body
    html += "<tbody>"
    t_row_group_names = list(group_to_sample_to_anchor_to_td.keys())
    if dt.pconfig.sort_rows:
        t_row_group_names = natsorted(t_row_group_names)

    non_trivial_groups_present = any(len(group_to_sample_to_anchor_to_td[g_name]) > 1 for g_name in t_row_group_names)

    for g_name in t_row_group_names:
        group_classes: List[str] = []
        # Hide the row if all cells are empty or hidden
        all_samples_empty = True
        for s_name in group_to_sample_to_anchor_to_td[g_name]:
            if not all(group_to_sample_to_anchor_to_empty[g_name][s_name].values()):  # not all empty!
                all_samples_empty = False
                break
        if all_samples_empty:
            group_classes.append("row-empty")
        for number_in_group, s_name in enumerate(group_to_sample_to_anchor_to_td[g_name]):
            tr_classes: List[str] = []
            prefix = ""
            if non_trivial_groups_present:
                caret_cls = ""
                if len(group_to_sample_to_anchor_to_td[g_name]) > 1 and number_in_group == 0:
                    caret_cls = "expandable-row-caret"
                    tr_classes.append("expandable-row-primary")
                prefix += f'<div style="display: inline-block; width: 20px" class="{caret_cls}">&nbsp;</div>'
            if number_in_group != 0:
                prefix += "&nbsp;↳&nbsp;"
                tr_classes.append("expandable-row-secondary expandable-row-secondary-hidden")
            cls = " ".join(group_classes + tr_classes)
            html += f'<tr data-sample-group="{escape(g_name)}" data-table-id="{dt.id}" class="{cls}">'
            # Sample name row header
            html += f'<th class="rowheader" data-sorting-val="{escape(g_name)}">{prefix}<span class="th-sample-name" data-original-sn="{escape(s_name)}">{s_name}</span></th>'
            for col_anchor in col_to_th.keys():
                cell_html = group_to_sample_to_anchor_to_td[g_name][s_name].get(col_anchor)
                if not cell_html:
                    td_hide_cls = "column-hidden" if col_to_hidden[col_anchor] else ""
                    sorting_val = group_to_sorting_to_anchor_to_val.get(g_name, {}).get(col_anchor, "")
                    cell_html = (
                        f'<td class="data-coloured {col_anchor} {td_hide_cls}" data-sorting-val="{sorting_val}"></td>'
                    )
                html += cell_html
            html += "</tr>"
    html += "</tbody></table></div>"
    if len(group_to_sample_to_anchor_to_td) > 10 and config.collapse_tables:
        html += (
            f'<div class="mqc-table-expand"><span>Expand table</span> {get_material_icon("mdi:chevron-down", 20)}</div>'
        )
    html += "</div>"

    # Save the raw values to a file if requested
    if dt.pconfig.save_file:
        fname = dt.pconfig.raw_data_fn or f"multiqc_{dt.anchor}"
        flatten_raw_vals: Dict[str, Dict[str, ValueT]] = {}
        for g_name, g_data in group_to_sample_to_anchor_to_val.items():
            for s_name, s_data in g_data.items():
                flatten_raw_vals[str(s_name)] = {str(k): v for k, v in s_data.items()}
        report.write_data_file(flatten_raw_vals, fname)
        if config.data_dump_file_write_raw:
            report.write_data_file(flatten_raw_vals, fname, data_format="json")
            report.saved_raw_data_keys[fname] = None

    # Build the bootstrap modal to customise columns and order
    modal = ""
    if not config.simple_output and add_control_panel and not _is_configure_columns_disabled(len(col_to_th)):
        modal = _configuration_modal(
            table_anchor=dt.anchor,
            title=table_title,
            trows="".join(col_to_modal_headers.values()),
            violin_anchor=violin_anchor,
        )

    return html, modal


def _configuration_modal(table_anchor: str, title: str, trows: str, violin_anchor: Optional[str] = None) -> str:
    data = f"data-table-anchor='{table_anchor}'"
    if violin_anchor is not None:
        data += f" data-violin-anchor='{violin_anchor}'"
    return f"""
    <!-- MultiQC Table Columns Modal -->
    <div class="modal fade mqc_config_modal" id="{table_anchor}_config_modal" tabindex="-1">
      <div class="modal-dialog modal-xl">
        <div class="modal-content">
          <div class="modal-header">
            <h4 class="modal-title">{title}: Columns</h4>
            <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="Close"></button>
          </div>
          <div class="modal-body">
            <p>Uncheck the tick box to hide columns. Click and drag the handle on the left to change order. Table ID: <code>{table_anchor}</code></p>
            <p>
                <button class="btn btn-outline-secondary btn-sm mqc_config_modal_bulk_visible" {data} data-action="showAll">Show All</button>
                <button class="btn btn-outline-secondary btn-sm mqc_config_modal_bulk_visible" {data} data-action="showNone">Show None</button>
            </p>
            <table class="table mqc_table mqc_sortable mqc_config_modal_table" id="{table_anchor}_config_modal_table" data-title="{title}">
              <thead>
                <tr>
                  <th class="sorthandle" style="text-align:center;">Sort</th>
                  <th style="text-align:center;">Visible</th>
                  <th>Group</th>
                  <th>Column</th>
                  <th>Description</th>
                  <th>ID</th>
                  <th>Scale</th>
                </tr>
              </thead>
              <tbody>
                {trows}
              </tbody>
            </table>
        </div>
        <div class="modal-footer"> <button type="button" class="btn btn-outline-secondary" data-bs-dismiss="modal">Close</button> </div>
    </div> </div> </div>"""


def _get_sortlist_js(dt: DataTable) -> str:
    """
    Custom column sorting order string for JavaScript table plot property.

    The order is provided in the following form:

    ```yaml
    custom_plot_config:
      general_stats_table:
        defaultsort:
          - column: "Mean Insert Length"
            direction: asc
          - column: "Starting Amount (ng)"
      quast_table:
        defaultsort:
        - column: "Largest contig"
    ```

    It is returned in a form os a list literal, as expected by the jQuery tablesorter plugin.
    """
    defaultsort = dt.pconfig.defaultsort
    if defaultsort is None:
        return ""

    headers = dt.get_headers_in_order()
    sortlist: List[Tuple[int, int]] = []

    # defaultsort is a list of {column, direction} objects
    for d in defaultsort:
        try:
            # The first element of the triple is not actually unique, it's a bucket index,
            # so we must re-enumerate ourselves here
            idx = next(
                idx
                for idx, (_, k, header) in enumerate(headers)
                if d["column"].lower() in [k.lower(), header.title.lower()]
            )
        except StopIteration:
            logger.warning(
                "Tried to sort by column '%s', but column was not found. Available columns: %s",
                d["column"],
                [k for (_, k, _) in headers],
            )
            return ""
        idx += 1  # to account for col1_header
        direction = 0 if d.get("direction", "").startswith("asc") else 1
        sortlist.append((idx, direction))

    sortlist_formatted = [list(t) for t in sortlist]
    return str(sortlist_formatted)


def _is_configure_columns_disabled(num_columns: int) -> bool:
    return num_columns > config.max_configurable_table_columns
