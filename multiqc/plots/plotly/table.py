import logging
from collections import defaultdict
from typing import Tuple, Optional, List

from multiqc.plots.table_object import DataTable
from multiqc.utils import config, mqc_colour, util_functions, report

logger = logging.getLogger(__name__)


def plot(dt: List[DataTable]) -> str:
    from multiqc.plots.plotly import violin

    return violin.plot(dt, show_table_by_default=True)


def make_table(dt: DataTable, violin_id: Optional[str] = None) -> Tuple[str, str]:
    """
    Build HTML for a MultiQC table, and HTML for the modal for configuring the table.
    :param dt: MultiQC datatable object
    :param violin_id: optional, will add a button to switch to a violin plot with this ID
    """

    t_headers = dict()
    t_modal_headers = dict()
    t_rows = dict()
    t_rows_empty = dict()
    dt.raw_vals = defaultdict(lambda: dict())
    empty_cells = dict()
    hidden_cols = 1
    table_title = dt.pconfig.get("table_title")
    if table_title is None:
        table_title = dt.id.replace("_", " ").title()

    for idx, k, header in dt.get_headers_in_order():
        rid = header["rid"]

        # Build the table header cell
        shared_key = ""
        if header.get("shared_key", None) is not None:
            shared_key = f" data-shared-key={header['shared_key']}"

        hide = ""
        muted = ""
        checked = ' checked="checked"'
        if header.get("hidden", False) is True:
            hide = "hidden"
            muted = " text-muted"
            checked = ""
            hidden_cols += 1

        data_attr = 'data-dmax="{}" data-dmin="{}" data-namespace="{}" {}'.format(
            header["dmax"], header["dmin"], header["namespace"], shared_key
        )

        ns = f'{header["namespace"]}: ' if header["namespace"] else ""
        cell_contents = f'<span class="mqc_table_tooltip" title="{ns}{header["description"]}">{header["title"]}</span>'

        t_headers[rid] = '<th id="header_{rid}" class="{rid} {h}" {da}>{c}</th>'.format(
            rid=rid, h=hide, da=data_attr, c=cell_contents
        )

        empty_cells[rid] = f'<td class="data-coloured {rid} {hide}"></td>'

        # Build the modal table row
        data = f"data-table-id='{dt.id}'"
        if violin_id:
            data += f" data-violin-id='{violin_id}'"
        t_modal_headers[rid] = f"""
        <tr class="{rid}{muted}" style="background-color: rgba({header["colour"]}, 0.15);">
          <td class="sorthandle ui-sortable-handle">||</span></td>
          <td style="text-align:center;">
            <input class="mqc_table_col_visible" type="checkbox" {checked} value="{rid}" {data}>
          </td>
          <td>{header["namespace"]}</td>
          <td>{header["title"]}</td>
          <td>{header["description"]}</td>
          <td><code>{k}</code></td>
          <td>{header.get("shared_key", "")}</td>
        </tr>"""

        # Make a colour scale
        if header["scale"] is False:
            c_scale = None
        else:
            c_scale = mqc_colour.mqc_colour_scale(
                name=header["scale"],
                minval=header["dmin"],
                maxval=header["dmax"],
                id=dt.id,
            )

        # Collect conditional formatting config
        cond_formatting_rules = {}
        if header.get("cond_formatting_rules"):
            cond_formatting_rules[rid] = header["cond_formatting_rules"]
        cond_formatting_rules.update(config.table_cond_formatting_rules)

        cond_formatting_colours = header.get("cond_formatting_colours", [])
        cond_formatting_colours.extend(config.table_cond_formatting_colours)

        # Add the data table cells
        for s_name, samp in dt.data[idx].items():
            if k in samp:
                val = samp[k]
                kname = f"{header['namespace']}_{rid}"
                dt.raw_vals[s_name][kname] = val

                if "modify" in header and callable(header["modify"]):
                    try:
                        val = header["modify"](val)
                    except TypeError as e:
                        logger.debug(f"Error modifying table value {kname} : {val} - {e}")

                if c_scale and c_scale.name not in c_scale.qualitative_scales:
                    try:
                        dmin = header["dmin"]
                        dmax = header["dmax"]
                        percentage = ((float(val) - dmin) / (dmax - dmin)) * 100
                        # Treat 0 as 0-width and make bars width of absolute value
                        if header.get("bars_zero_centrepoint"):
                            dmax = max(abs(header["dmin"]), abs(header["dmax"]))
                            dmin = 0
                            percentage = ((abs(float(val)) - dmin) / (dmax - dmin)) * 100
                        percentage = min(percentage, 100)
                        percentage = max(percentage, 0)
                    except (ZeroDivisionError, ValueError, TypeError):
                        percentage = 0
                else:
                    percentage = 100

                if "format" in header and callable(header["format"]):
                    valstring = header["format"](val)
                else:
                    try:
                        # "format" is a format string?
                        valstring = str(header["format"].format(val))
                    except ValueError:
                        try:
                            valstring = str(header["format"].format(float(val)))
                        except ValueError:
                            valstring = str(val)
                    except Exception:
                        valstring = str(val)

                    # This is horrible, but Python locale settings are worse
                    if config.thousandsSep_format is None:
                        config.thousandsSep_format = '<span class="mqc_small_space"></span>'
                    if config.decimalPoint_format is None:
                        config.decimalPoint_format = "."
                    valstring = valstring.replace(".", "DECIMAL").replace(",", "THOUSAND")
                    valstring = valstring.replace("DECIMAL", config.decimalPoint_format).replace(
                        "THOUSAND", config.thousandsSep_format
                    )

                suffix = header.get("suffix")
                if suffix:
                    # Add a space before the suffix, but not as an actual character, so ClipboardJS would copy
                    # the whole value without the space. Also, remove &nbsp; that we don't want ClipboardJS to copy.
                    suffix = suffix.replace("&nbsp;", " ").strip()
                    valstring += "<span class='mqc_small_space'></span>" + suffix

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
                                    if "eq" in cmp and float(cmp["eq"]) == float(val):
                                        cmatches[ftype] = True
                                    if "ne" in cmp and float(cmp["ne"]) != float(val):
                                        cmatches[ftype] = True
                                    if "gt" in cmp and float(cmp["gt"]) < float(val):
                                        cmatches[ftype] = True
                                    if "lt" in cmp and float(cmp["lt"]) > float(val):
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
                    valstring = f'<span class="badge" style="background-color:{badge_col}">{valstring}</span>'

                # Determine background color based on scale. Only relevant for hashable values. If value is for some
                # reason a dict or a list, it's not hashable and the logic determining the color will not work.
                hashable = True
                try:
                    hash(val)
                except TypeError:
                    hashable = False
                    print(f"Value {val} is not hashable for table {dt.id}, column {k}, sample {s_name}")

                # Categorical background colours supplied
                if hashable and val in header.get("bgcols", {}).keys():
                    col = f"style=\"background-color:{header['bgcols'][val]} !important;\""
                    if s_name not in t_rows:
                        t_rows[s_name] = dict()
                    t_rows[s_name][rid] = f'<td val="{val}" class="{rid} {hide}" {col}>{valstring}</td>'

                # Build table cell background colour bar
                elif hashable and header["scale"]:
                    if c_scale is not None:
                        col = " background-color:{} !important;".format(
                            c_scale.get_colour(val, source=f'Table "{dt.id}", column "{k}"')
                        )
                    else:
                        col = ""
                    bar_html = f'<span class="bar" style="width:{percentage}%;{col}"></span>'
                    val_html = f'<span class="val">{valstring}</span>'
                    wrapper_html = f'<div class="wrapper">{bar_html}{val_html}</div>'

                    if s_name not in t_rows:
                        t_rows[s_name] = dict()
                    t_rows[s_name][rid] = f'<td val="{val}" class="data-coloured {rid} {hide}">{wrapper_html}</td>'

                # Scale / background colours are disabled
                else:
                    if s_name not in t_rows:
                        t_rows[s_name] = dict()
                    t_rows[s_name][rid] = f'<td val="{val}" class="{rid} {hide}">{valstring}</td>'

                # Is this cell hidden or empty?
                if s_name not in t_rows_empty:
                    t_rows_empty[s_name] = dict()
                t_rows_empty[s_name][rid] = header.get("hidden", False) or str(val).strip() == ""

        # Remove header if we don't have any filled cells for it
        if sum([len(rows) for rows in t_rows.values()]) == 0:
            if header.get("hidden", False) is True:
                hidden_cols -= 1
            t_headers.pop(rid, None)
            t_modal_headers.pop(rid, None)
            logger.debug(f"Removing header {k} from table, as no data")

    #
    # Put everything together
    #

    # Buttons above the table
    html = ""
    if not config.simple_output:
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
            buttons.append(
                f"""
            <button type="button" class="mqc_table_configModal_btn btn btn-default btn-sm" data-toggle="modal" 
                data-target="#{dt.id}_configModal">
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
    col1_header = dt.pconfig.get("col1_header", "Sample Name")
    html += f"<thead><tr><th class=\"rowheader\">{col1_header}</th>{''.join(t_headers.values())}</tr></thead>"

    # Build the table body
    html += "<tbody>"
    t_row_keys = t_rows.keys()
    if dt.pconfig.get("sort_rows", dt.pconfig.get("sortRows")) is not False:
        t_row_keys = sorted(t_row_keys)
    for s_name in t_row_keys:
        # Hide the row if all cells are empty or hidden
        row_hidden = ' style="display:none"' if all(t_rows_empty[s_name].values()) else ""
        html += f"<tr{row_hidden}>"
        # Sample name row header
        html += f'<th class="rowheader" data-original-sn="{s_name}">{s_name}</th>'
        for k in t_headers:
            html += t_rows[s_name].get(k, empty_cells[k])
        html += "</tr>"
    html += "</tbody></table></div>"
    if len(t_rows) > 10 and config.collapse_tables:
        html += '<div class="mqc-table-expand"><span class="glyphicon glyphicon-chevron-down" aria-hidden="true"></span></div>'
    html += "</div>"

    # Save the raw values to a file if requested
    if dt.pconfig.get("save_file") is True:
        fn = dt.pconfig.get("raw_data_fn", f"multiqc_{dt.id}")
        util_functions.write_data_file(dt.raw_vals, fn)
        report.saved_raw_data[fn] = dt.raw_vals

    # Build the bootstrap modal to customise columns and order
    modal = ""
    if not config.simple_output:
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
    <div class="modal fade" id="{tid}_configModal" tabindex="-1">
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
    defaultsort = dt.pconfig.get("defaultsort")
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
                if d["column"].lower() in [k.lower(), header["title"].lower()]
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
