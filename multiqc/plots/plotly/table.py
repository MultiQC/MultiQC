import logging
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Union

from multiqc import config, report
from multiqc.plots.table_object import ColumnAnchor, DataTable, SampleGroup, SampleName, ValueT
from multiqc.utils import mqc_colour
from typing import TYPE_CHECKING
from natsort import natsorted

if TYPE_CHECKING:  # to avoid circular import
    from multiqc.plots.plotly.violin import ViolinPlot

logger = logging.getLogger(__name__)


def plot(dt: List[DataTable]) -> "ViolinPlot":
    from multiqc.plots.plotly import violin

    return violin.plot(dt, show_table_by_default=True)


def make_table(
    dt: DataTable,
    violin_anchor: Optional[str] = None,
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
        col_anchor: ColumnAnchor = header.rid

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
        section = dt.sections[idx]
        for group_name, group_rows in section.rows_by_sgroup.items():
            for row_idx, row in enumerate(group_rows):
                if col_key not in row.raw_data:
                    continue

                val: ValueT = row.raw_data[col_key]
                valstr: str = row.formatted_data[col_key]

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

    #
    # Put everything together
    #

    html = ""

    # Buttons above the table
    if not config.simple_output and add_control_panel:
        # Copy Table Button
        buttons: List[str] = []

        buttons.append(
            f"""
        <button type="button" class="mqc_table_copy_btn btn btn-default btn-sm" data-clipboard-target="table#{dt.anchor}">
            <span class="glyphicon glyphicon-copy"></span> Copy table
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
            <button type="button" class="mqc_table_config_modal_btn btn btn-default btn-sm {disabled_class}" data-toggle="modal"
                data-target="#{dt.anchor}_config_modal" {disabled_attrs}>
                <span class="glyphicon glyphicon-th"></span> Configure columns
            </button>
            """
            )

        # Sort By Highlight button
        buttons.append(
            f"""
        <button type="button" class="mqc_table_sortHighlight btn btn-default btn-sm"
            data-table-anchor="{dt.anchor}" data-direction="desc" style="display:none;">
            <span class="glyphicon glyphicon-sort-by-attributes-alt"></span> Sort by highlight
        </button>
        """
        )

        # Scatter Plot Button
        if len(col_to_th) > 1:
            buttons.append(
                f"""
            <button type="button" class="mqc_table_make_scatter btn btn-default btn-sm"
                data-toggle="modal" data-target="#table_scatter_modal" data-table-anchor="{dt.anchor}">
                <span class="glyphicon glyphicon glyphicon-equalizer"></span> Scatter plot
            </button>
            """
            )

        if violin_anchor is not None:
            buttons.append(
                f"""
            <button type="button" class="mqc-table-to-violin btn btn-default btn-sm"
                data-table-anchor="{dt.anchor}" data-violin-anchor="{violin_anchor}">
                <span class="glyphicon glyphicon-align-left"></span> Violin plot
            </button>
            """
            )

        buttons.append(
            f"""
        <button type="button" class="export-plot btn btn-default btn-sm"
            data-plot-anchor="{violin_anchor or dt.anchor}" data-type="table"
        >Export as CSV</button>
        """
        )

        # "Showing x of y columns" text
        row_visibilities = [
            all(group_to_sample_to_anchor_to_empty[s_name].values()) for s_name in group_to_sample_to_anchor_to_empty
        ]
        visible_rows = [x for x in row_visibilities if not x]

        # Visible rows
        t_showing_rows_txt = f'Showing <sup id="{dt.anchor}_numrows" class="mqc_table_numrows">{len(visible_rows)}</sup>/<sub>{len(group_to_sample_to_anchor_to_td)}</sub> rows'

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

        panel = "\n".join(buttons)
        html += f"""
        <div class='row'>\n<div class='col-xs-12'>\n{panel}\n</div>\n</div>
        """

    # Build the table itself
    collapse_class = (
        "mqc-table-collapse" if len(group_to_sample_to_anchor_to_td) > 10 and config.collapse_tables else ""
    )
    html += f"""
        <div id="{dt.anchor}_container" class="mqc_table_container">
            <div class="table-responsive mqc-table-responsive {collapse_class}">
                <table id="{dt.anchor}" class="table table-condensed mqc_table mqc_per_sample_table" data-title="{table_title}" data-sortlist="{_get_sortlist(dt)}">
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
        html += '<div class="mqc-table-expand"><span class="glyphicon glyphicon-chevron-down" aria-hidden="true"></span></div>'
    html += "</div>"

    # Save the raw values to a file if requested
    if dt.pconfig.save_file:
        fname = dt.pconfig.raw_data_fn or f"multiqc_{dt.anchor}"
        flatten_raw_vals: Dict[str, Dict[str, ValueT]] = {}
        for g_name, g_data in group_to_sample_to_anchor_to_val.items():
            for s_name, s_data in g_data.items():
                flatten_raw_vals[str(s_name)] = {str(k): v for k, v in s_data.items()}
        report.write_data_file(flatten_raw_vals, fname)
        report.saved_raw_data[fname] = flatten_raw_vals

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
      <div class="modal-dialog modal-lg">
        <div class="modal-content">
          <div class="modal-header">
            <button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
            <h4 class="modal-title">{title}: Columns</h4>
          </div>
          <div class="modal-body">
            <p>Uncheck the tick box to hide columns. Click and drag the handle on the left to change order. Table ID: <code>{table_anchor}</code></p>
            <p>
                <button class="btn btn-default btn-sm mqc_config_modal_bulk_visible" {data} data-action="showAll">Show All</button>
                <button class="btn btn-default btn-sm mqc_config_modal_bulk_visible" {data} data-action="showNone">Show None</button>
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
        <div class="modal-footer"> <button type="button" class="btn btn-default" data-dismiss="modal">Close</button> </div>
    </div> </div> </div>"""


def _get_sortlist(dt: DataTable) -> str:
    """
    Custom column sorting order for a table plot. The order is provided in the following form:

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

    return str(sortlist)


def _is_configure_columns_disabled(num_columns: int) -> bool:
    return num_columns > 50
