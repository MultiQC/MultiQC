import logging
from collections import defaultdict
from typing import Tuple, Optional, List, Dict

from multiqc.plots.table_object import DataTable, ValueT
from multiqc import config, report
from multiqc.utils import mqc_colour

logger = logging.getLogger(__name__)


def plot(dt: List[DataTable]):
    from multiqc.plots.plotly import violin

    return violin.plot(dt, show_table_by_default=True)


def make_table(
    dt: DataTable,
    violin_id: Optional[str] = None,
    add_control_panel: bool = True,
) -> Tuple[str, str]:
    """
    Build HTML for a MultiQC table, and HTML for the modal for configuring the table.
    :param dt: MultiQC datatable object
    :param violin_id: optional, will add a button to switch to a violin plot with this ID
    :param add_control_panel: whether to add the control panel with buttons above the table
    """

    t_headers = dict()
    t_modal_headers = dict()
    t_rows: Dict[str, Dict[str, str]] = dict()
    t_rows_empty: Dict[str, Dict[str, bool]] = dict()
    raw_vals: Dict[str, Dict[str, ValueT]] = defaultdict(lambda: dict())
    empty_cells = dict()
    hidden_cols = 1
    table_title = dt.pconfig.title
    if table_title is None:
        table_title = dt.id.replace("_", " ").title()

    def escape(s: str) -> str:
        return s.replace('"', "&quot;").replace("'", "&#39;").replace("<", "&lt;").replace(">", "&gt;")

    for idx, metric, header in dt.get_headers_in_order():
        rid = header.rid

        # Build the table header cell
        shared_key = ""
        if header.shared_key is not None:
            shared_key = f" data-shared-key={header.shared_key}"

        hide = ""
        muted = ""
        checked = ' checked="checked"'
        if header.hidden:
            hide = "hidden"
            muted = " text-muted"
            checked = ""
            hidden_cols += 1

        data_attr = 'data-dmax="{}" data-dmin="{}" data-namespace="{}" {}'.format(
            header.dmax, header.dmin, header.namespace, shared_key
        )

        ns = f"{header.namespace}: " if header.namespace else ""
        cell_contents = (
            f'<span class="mqc_table_tooltip" title="{ns}{header.description}" data-html="true">{header.title}</span>'
        )

        t_headers[rid] = '<th id="header_{rid}" class="{rid} {h}" {da}>{c}</th>'.format(
            rid=rid, h=hide, da=data_attr, c=cell_contents
        )

        empty_cells[rid] = f'<td class="data-coloured {rid} {hide}"></td>'

        # Build the modal table row
        data = f"data-table-id='{dt.id}'"
        if violin_id:
            data += f" data-violin-id='{violin_id}'"
        t_modal_headers[rid] = f"""
        <tr class="{rid}{muted}" style="background-color: rgba({header.color}, 0.15);">
          <td class="sorthandle ui-sortable-handle">||</span></td>
          <td style="text-align:center;">
            <input class="mqc_table_col_visible" type="checkbox" {checked} value="{rid}" {data}>
          </td>
          <td>{header.namespace}</td>
          <td>{header.title}</td>
          <td>{header.description}</td>
          <td><code>{metric}</code></td>
          <td>{header.shared_key or ""}</td>
        </tr>"""

        # Make a colour scale
        if header.scale is False:
            c_scale = None
        else:
            c_scale = mqc_colour.mqc_colour_scale(
                name=header.scale,
                minval=header.dmin,
                maxval=header.dmax,
                id=dt.id,
            )

        # Collect conditional formatting config
        cond_formatting_rules = {}
        if header.cond_formatting_rules:
            cond_formatting_rules[rid] = header.cond_formatting_rules
        cond_formatting_rules.update(config.table_cond_formatting_rules)

        cond_formatting_colours = header.cond_formatting_colours
        cond_formatting_colours.extend(config.table_cond_formatting_colours)

        # Add the data table cells
        for s_name in dt.raw_data[idx].keys():
            if metric in dt.raw_data[idx][s_name]:
                val: ValueT = dt.raw_data[idx][s_name][metric]
                valstr: str = dt.formatted_data[idx][s_name][metric]

                raw_vals[s_name][f"{header.namespace}_{rid}"] = val

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
                valstr = valstr

                suffix = header.suffix
                if suffix:
                    # Add a space before the suffix, but not as an actual character, so ClipboardJS would copy
                    # the whole value without the space. Also, remove &nbsp; that we don't want ClipboardJS to copy.
                    suffix = suffix.replace("&nbsp;", " ").strip()
                    valstr += "<span class='mqc_small_space'></span>" + suffix

                # Conditional formatting
                # Build empty dict for cformatting matches
                cmatches = {}
                for cfc in cond_formatting_colours:
                    for cfck in cfc:
                        cmatches[cfck] = False
                # Find general rules followed by column-specific rules
                for cfk in ["all_columns", rid, dt.id]:
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
                                except Exception:
                                    logger.warning(f"Not able to apply table conditional formatting to '{val}' ({cmp})")
                # Apply HTML in order of config keys
                badge_col = None
                for cfc in cond_formatting_colours:
                    for cfck in cfc:  # should always be one, but you never know
                        if cmatches[cfck]:
                            badge_col = cfc[cfck]
                if badge_col is not None:
                    valstr = f'<span class="badge" style="background-color:{badge_col}">{valstr}</span>'

                # Determine background color based on scale. Only relevant for hashable values. If value is for some
                # reason a dict or a list, it's not hashable and the logic determining the color will not work.
                hashable = True
                try:
                    hash(val)
                except TypeError:
                    hashable = False
                    print(f"Value {val} is not hashable for table {dt.id}, column {metric}, sample {s_name}")

                # Categorical background colours supplied
                if isinstance(val, str) and val in header.bgcols.keys():
                    col = f'style="background-color:{header.bgcols[val]} !important;"'
                    if s_name not in t_rows:
                        t_rows[s_name] = dict()
                    t_rows[s_name][rid] = f'<td val="{escape(str(val))}" class="{rid} {hide}" {col}>{valstr}</td>'

                # Build table cell background colour bar
                elif hashable and header.scale:
                    if c_scale is not None:
                        col = " background-color:{} !important;".format(
                            c_scale.get_colour(val, source=f'Table "{dt.id}", column "{metric}"')
                        )
                    else:
                        col = ""
                    bar_html = f'<span class="bar" style="width:{percentage}%;{col}"></span>'
                    val_html = f'<span class="val">{valstr}</span>'
                    wrapper_html = f'<div class="wrapper">{bar_html}{val_html}</div>'

                    if s_name not in t_rows:
                        t_rows[s_name] = dict()
                    t_rows[s_name][rid] = (
                        f'<td val="{escape(str(val))}" class="data-coloured {rid} {hide}">{wrapper_html}</td>'
                    )

                # Scale / background colours are disabled
                else:
                    if s_name not in t_rows:
                        t_rows[s_name] = dict()
                    t_rows[s_name][rid] = f'<td val="{escape(str(val))}" class="{rid} {hide}">{valstr}</td>'

                # Is this cell hidden or empty?
                if s_name not in t_rows_empty:
                    t_rows_empty[s_name] = dict()
                t_rows_empty[s_name][rid] = header.hidden or str(val).strip() == ""

        # Remove header if we don't have any filled cells for it
        if sum([len(rows) for rows in t_rows.values()]) == 0:
            if header.hidden:
                hidden_cols -= 1
            t_headers.pop(rid, None)
            t_modal_headers.pop(rid, None)
            logger.debug(f"Removing header {metric} from table, as no data")

    #
    # Put everything together
    #

    # Buttons above the table
    html = ""
    if not config.simple_output and add_control_panel:
        # Copy Table Button
        buttons = []

        buttons.append(
            f"""
        <button type="button" class="mqc_table_copy_btn btn btn-default btn-sm" data-clipboard-target="table#{dt.id}">
            <span class="glyphicon glyphicon-copy"></span> Copy table
        </button>
        """
        )

        # Configure Columns Button
        if len(t_headers) > 1:
            # performance degrades substantially when configuring thousands of columns
            # it is effectively unusable.
            disabled_class = ""
            disabled_attrs = ""
            if _is_configure_columns_disabled(len(t_headers)):
                disabled_class = "mqc_table_tooltip"
                disabled_attrs = 'disabled title="Table is too large to configure columns"'

            buttons.append(
                f"""
            <button type="button" class="mqc_table_configModal_btn btn btn-default btn-sm {disabled_class}" data-toggle="modal"
                data-target="#{dt.id}_configModal" {disabled_attrs}>
                <span class="glyphicon glyphicon-th"></span> Configure columns
            </button>
            """
            )

        # Sort By Highlight button
        buttons.append(
            f"""
        <button type="button" class="mqc_table_sortHighlight btn btn-default btn-sm"
            data-target="#{dt.id}" data-direction="desc" style="display:none;">
            <span class="glyphicon glyphicon-sort-by-attributes-alt"></span> Sort by highlight
        </button>
        """
        )

        # Scatter Plot Button
        if len(t_headers) > 1:
            buttons.append(
                f"""
            <button type="button" class="mqc_table_makeScatter btn btn-default btn-sm"
                data-toggle="modal" data-target="#tableScatterModal" data-table="#{dt.id}">
                <span class="glyphicon glyphicon glyphicon-equalizer"></span> Scatter plot
            </button>
            """
            )

        if violin_id is not None:
            buttons.append(
                f"""
            <button type="button" class="mqc-table-to-violin btn btn-default btn-sm"
                data-table-id="{dt.id}" data-violin-id="{violin_id}">
                <span class="glyphicon glyphicon-align-left"></span> Violin plot
            </button>
            """
            )

        buttons.append(
            f"""
        <button type="button" class="export-plot btn btn-default btn-sm"
            data-pid="{violin_id or dt.id}" data-type="table"
        >Export as CSV</button>
        """
        )

        # "Showing x of y columns" text
        row_visibilities = [all(t_rows_empty[s_name].values()) for s_name in t_rows_empty]
        visible_rows = [x for x in row_visibilities if not x]

        # Visible rows
        t_showing_rows_txt = f'Showing <sup id="{dt.id}_numrows" class="mqc_table_numrows">{len(visible_rows)}</sup>/<sub>{len(t_rows)}</sub> rows'

        # How many columns are visible?
        ncols_vis = (len(t_headers) + 1) - hidden_cols
        t_showing_cols_txt = ""
        if len(t_headers) > 1:
            t_showing_cols_txt = f' and <sup id="{dt.id}_numcols" class="mqc_table_numcols">{ncols_vis}</sup>/<sub>{len(t_headers)}</sub> columns'

        # Build table header text
        buttons.append(
            f"""
        <small id="{dt.id}_numrows_text" class="mqc_table_numrows_text">{t_showing_rows_txt}{t_showing_cols_txt}.</small>
        """
        )

        panel = "\n".join(buttons)
        html += f"""
        <div class='row'>\n<div class='col-xs-12'>\n{panel}\n</div>\n</div>
        """

    # Build the table itself
    collapse_class = "mqc-table-collapse" if len(t_rows) > 10 and config.collapse_tables else ""
    html += f"""
        <div id="{dt.id}_container" class="mqc_table_container">
            <div class="table-responsive mqc-table-responsive {collapse_class}">
                <table id="{dt.id}" class="table table-condensed mqc_table" data-title="{table_title}" data-sortlist="{_get_sortlist(dt)}">
        """

    # Build the header row
    col1_header = dt.pconfig.col1_header
    html += f"<thead><tr><th class=\"rowheader\">{col1_header}</th>{''.join(t_headers.values())}</tr></thead>"

    # Build the table body
    html += "<tbody>"
    t_row_keys = list(t_rows.keys())
    if dt.pconfig.sort_rows:
        t_row_keys = sorted(t_row_keys)
    for s_name in t_row_keys:
        # Hide the row if all cells are empty or hidden
        row_hidden = ' style="display:none"' if all(t_rows_empty[s_name].values()) else ""
        html += f"<tr{row_hidden}>"
        # Sample name row header
        html += f'<th class="rowheader" data-original-sn="{escape(s_name)}">{s_name}</th>'
        for metric in t_headers:
            html += t_rows[s_name].get(metric, empty_cells[metric])
        html += "</tr>"
    html += "</tbody></table></div>"
    if len(t_rows) > 10 and config.collapse_tables:
        html += '<div class="mqc-table-expand"><span class="glyphicon glyphicon-chevron-down" aria-hidden="true"></span></div>'
    html += "</div>"

    # Save the raw values to a file if requested
    if dt.pconfig.save_file:
        fn = dt.pconfig.raw_data_fn or f"multiqc_{dt.id}"
        report.write_data_file(raw_vals, fn)
        report.saved_raw_data[fn] = raw_vals

    # Build the bootstrap modal to customise columns and order
    modal = ""
    if not config.simple_output and add_control_panel and not _is_configure_columns_disabled(len(t_headers)):
        modal = _configuration_modal(
            tid=dt.id,
            title=table_title,
            trows="".join(t_modal_headers.values()),
            violin_id=violin_id,
        )

    return html, modal


def _configuration_modal(tid: str, title: str, trows: str, violin_id: Optional[str] = None) -> str:
    data = f"data-table-id='{tid}'"
    if violin_id is not None:
        data += f" data-violin-id='{violin_id}'"
    return f"""
    <!-- MultiQC Table Columns Modal -->
    <div class="modal fade mqc_configModal" id="{tid}_configModal" tabindex="-1">
      <div class="modal-dialog modal-lg">
        <div class="modal-content">
          <div class="modal-header">
            <button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
            <h4 class="modal-title">{title}: Columns</h4>
          </div>
          <div class="modal-body">
            <p>Uncheck the tick box to hide columns. Click and drag the handle on the left to change order. Table ID: <code>{tid}</code></p>
            <p>
                <button class="btn btn-default btn-sm mqc_configModal_bulkVisible" {data} data-action="showAll">Show All</button>
                <button class="btn btn-default btn-sm mqc_configModal_bulkVisible" {data} data-action="showNone">Show None</button>
            </p>
            <table class="table mqc_table mqc_sortable mqc_configModal_table" id="{tid}_configModal_table" data-title="{title}">
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
    sortlist = []

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
        sortlist.append([idx, direction])

    return str(sortlist)


def _is_configure_columns_disabled(num_columns: int) -> bool:
    return num_columns > 50
