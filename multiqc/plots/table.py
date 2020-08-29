#!/usr/bin/env python

""" MultiQC functions to plot a table """

from collections import defaultdict, OrderedDict
import logging
import random

from multiqc.utils import config, report, util_functions, mqc_colour
from multiqc.plots import table_object, beeswarm
logger = logging.getLogger(__name__)

letters = 'abcdefghijklmnopqrstuvwxyz'

def plot (data, headers=None, pconfig=None):
    """ Return HTML for a MultiQC table.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param headers: list of optional dicts with column config in key:value pairs.
    :return: HTML ready to be inserted into the page
    """
    if headers is None:
        headers = []
    if pconfig is None:
        pconfig = {}

    # Allow user to overwrite any given config for this plot
    if 'id' in pconfig and pconfig['id'] and pconfig['id'] in config.custom_plot_config:
        for k, v in config.custom_plot_config[pconfig['id']].items():
            pconfig[k] = v

    # Make a datatable object
    dt = table_object.datatable(data, headers, pconfig)

    # Collect unique sample names
    s_names = set()
    for d in dt.data:
        for s_name in d.keys():
            s_names.add(s_name)

    # Make a beeswarm plot if we have lots of samples
    if len(s_names) >= config.max_table_rows and pconfig.get('no_beeswarm') is not True:
        logger.debug('Plotting beeswarm instead of table, {} samples'.format(len(s_names)))
        warning = '<p class="text-muted"><span class="glyphicon glyphicon-exclamation-sign" ' \
            'title="A beeswarm plot has been generated instead because of the large number of samples. '\
            'See http://multiqc.info/docs/#tables--beeswarm-plots"'\
            ' data-toggle="tooltip"></span> Showing {} samples.</p>'.format(len(s_names))
        return warning + beeswarm.make_plot( dt )
    else:
        return make_table ( dt )


def make_table (dt):
    """
    Build the HTML needed for a MultiQC table.
    :param data: MultiQC datatable object
    """

    table_id = dt.pconfig.get('id', 'table_{}'.format(''.join(random.sample(letters, 4))) )
    table_id = report.save_htmlid(table_id)
    t_headers = OrderedDict()
    t_modal_headers = OrderedDict()
    t_rows = OrderedDict()
    t_rows_empty = OrderedDict()
    dt.raw_vals = defaultdict(lambda: dict())
    empty_cells = dict()
    hidden_cols = 1
    table_title = dt.pconfig.get('table_title')
    if table_title is None:
        table_title = table_id.replace("_", " ").title()

    for idx, k, header in dt.get_headers_in_order():

        rid = header['rid']

        # Build the table header cell
        shared_key = ''
        if header.get('shared_key', None) is not None:
            shared_key = ' data-shared-key={}'.format(header['shared_key'])

        hide = ''
        muted = ''
        checked = ' checked="checked"'
        if header.get('hidden', False) is True:
            hide = 'hidden'
            muted = ' text-muted'
            checked = ''
            hidden_cols += 1

        data_attr = 'data-dmax="{}" data-dmin="{}" data-namespace="{}" {}' \
            .format(header['dmax'], header['dmin'], header['namespace'], shared_key)

        cell_contents = '<span class="mqc_table_tooltip" title="{}: {}">{}</span>' \
            .format(header['namespace'], header['description'], header['title'])

        t_headers[rid] = '<th id="header_{rid}" class="{rid} {h}" {da}>{c}</th>' \
            .format(rid=rid, h=hide, da=data_attr, c=cell_contents)

        empty_cells[rid] = '<td class="data-coloured {rid} {h}"></td>'.format(rid=rid, h=hide)

        # Build the modal table row
        t_modal_headers[rid] = """
        <tr class="{rid}{muted}" style="background-color: rgba({col}, 0.15);">
          <td class="sorthandle ui-sortable-handle">||</span></td>
          <td style="text-align:center;">
            <input class="mqc_table_col_visible" type="checkbox" {checked} value="{rid}" data-target="#{tid}">
          </td>
          <td>{name}</td>
          <td>{title}</td>
          <td>{desc}</td>
          <td>{col_id}</td>
          <td>{sk}</td>
        </tr>""".format(
                rid = rid,
                muted = muted,
                checked = checked,
                tid = table_id,
                col = header['colour'],
                name = header['namespace'],
                title = header['title'],
                desc = header['description'],
                col_id = '<code>{}</code>'.format(k),
                sk = header.get('shared_key', '')
            )

        # Make a colour scale
        if header['scale'] == False:
            c_scale = None
        else:
            c_scale = mqc_colour.mqc_colour_scale(header['scale'], header['dmin'], header['dmax'])

        # Add the data table cells
        for (s_name, samp) in dt.data[idx].items():
            if k in samp:
                val = samp[k]
                kname = '{}_{}'.format(header['namespace'], rid)
                dt.raw_vals[s_name][kname] = val

                if 'modify' in header and callable(header['modify']):
                    val = header['modify'](val)

                try:
                    dmin = header['dmin']
                    dmax = header['dmax']
                    percentage = ((float(val) - dmin) / (dmax - dmin)) * 100
                    percentage = min(percentage, 100)
                    percentage = max(percentage, 0)
                except (ZeroDivisionError,ValueError):
                    percentage = 0

                try:
                    valstring = str(header['format'].format(val))
                except ValueError:
                    try:
                        valstring = str(header['format'].format(float(val)))
                    except ValueError:
                        valstring = str(val)
                except:
                    valstring = str(val)

                # This is horrible, but Python locale settings are worse
                if config.thousandsSep_format is None:
                    config.thousandsSep_format = '<span class="mqc_thousandSep"></span>'
                if config.decimalPoint_format is None:
                    config.decimalPoint_format = '.'
                valstring = valstring.replace('.', 'DECIMAL').replace(',', 'THOUSAND')
                valstring = valstring.replace('DECIMAL', config.decimalPoint_format).replace('THOUSAND', config.thousandsSep_format)

                # Percentage suffixes etc
                valstring += header.get('suffix', '')

                # Conditional formatting
                cmatches = { cfck: False for cfc in config.table_cond_formatting_colours for cfck in cfc }
                # Find general rules followed by column-specific rules
                for cfk in ['all_columns', rid]:
                    if cfk in config.table_cond_formatting_rules:
                        # Loop through match types
                        for ftype in cmatches.keys():
                            # Loop through array of comparison types
                            for cmp in config.table_cond_formatting_rules[cfk].get(ftype, []):
                                try:
                                    # Each comparison should be a dict with single key: val
                                    if 's_eq' in cmp and str(cmp['s_eq']).lower() == str(val).lower():
                                        cmatches[ftype] = True
                                    if 's_contains' in cmp and str(cmp['s_contains']).lower() in str(val).lower():
                                        cmatches[ftype] = True
                                    if 's_ne' in cmp and str(cmp['s_ne']).lower() != str(val).lower():
                                        cmatches[ftype] = True
                                    if 'eq' in cmp and float(cmp['eq']) == float(val):
                                        cmatches[ftype] = True
                                    if 'ne' in cmp and float(cmp['ne']) != float(val):
                                        cmatches[ftype] = True
                                    if 'gt' in cmp and float(cmp['gt']) < float(val):
                                        cmatches[ftype] = True
                                    if 'lt' in cmp and float(cmp['lt']) > float(val):
                                        cmatches[ftype] = True
                                except:
                                    logger.warning("Not able to apply table conditional formatting to '{}' ({})".format(val, cmp))
                # Apply HTML in order of config keys
                bgcol = None
                for cfc in config.table_cond_formatting_colours:
                    for cfck in cfc: # should always be one, but you never know
                        if cmatches[cfck]:
                            bgcol = cfc[cfck]
                if bgcol is not None:
                    valstring = '<span class="badge" style="background-color:{}">{}</span>'.format(bgcol, valstring)

                # Build HTML
                if not header['scale']:
                    if s_name not in t_rows:
                        t_rows[s_name] = dict()
                    t_rows[s_name][rid] = '<td class="{rid} {h}">{v}</td>'.format(rid=rid, h=hide, v=valstring)
                else:
                    if c_scale is not None:
                        col = ' background-color:{};'.format(c_scale.get_colour(val))
                    else:
                        col = ''
                    bar_html = '<span class="bar" style="width:{}%;{}"></span>'.format(percentage, col)
                    val_html = '<span class="val">{}</span>'.format(valstring)
                    wrapper_html = '<div class="wrapper">{}{}</div>'.format(bar_html, val_html)

                    if s_name not in t_rows:
                        t_rows[s_name] = dict()
                    t_rows[s_name][rid] = '<td class="data-coloured {rid} {h}">{c}</td>'.format(rid=rid, h=hide, c=wrapper_html)

                # Is this cell hidden or empty?
                if s_name not in t_rows_empty:
                    t_rows_empty[s_name] = dict()
                t_rows_empty[s_name][rid] = header.get('hidden', False) or str(val).strip() == ''

        # Remove header if we don't have any filled cells for it
        if sum([len(rows) for rows in t_rows.values()]) == 0:
            if header.get('hidden', False) is True:
                hidden_cols -= 1
            t_headers.pop(rid, None)
            t_modal_headers.pop(rid, None)
            logger.debug('Removing header {} from table, as no data'.format(k))

    #
    # Put everything together
    #

    # Buttons above the table
    html = ''
    if not config.simple_output:

        # Copy Table Button
        html += """
        <button type="button" class="mqc_table_copy_btn btn btn-default btn-sm" data-clipboard-target="#{tid}">
            <span class="glyphicon glyphicon-copy"></span> Copy table
        </button>
        """.format(tid=table_id)

        # Configure Columns Button
        if len(t_headers) > 1:
            html += """
            <button type="button" class="mqc_table_configModal_btn btn btn-default btn-sm" data-toggle="modal" data-target="#{tid}_configModal">
                <span class="glyphicon glyphicon-th"></span> Configure Columns
            </button>
            """.format(tid=table_id)

        # Sort By Highlight button
        html += """
        <button type="button" class="mqc_table_sortHighlight btn btn-default btn-sm" data-target="#{tid}" data-direction="desc" style="display:none;">
            <span class="glyphicon glyphicon-sort-by-attributes-alt"></span> Sort by highlight
        </button>
        """.format(tid=table_id)

        # Scatter Plot Button
        if len(t_headers) > 1:
            html += """
            <button type="button" class="mqc_table_makeScatter btn btn-default btn-sm" data-toggle="modal" data-target="#tableScatterModal" data-table="#{tid}">
                <span class="glyphicon glyphicon glyphicon-stats"></span> Plot
            </button>
            """.format(tid=table_id)

        # "Showing x of y columns" text
        row_visibilities = [ all(t_rows_empty[s_name].values()) for s_name in t_rows_empty ]
        visible_rows = [ x for x in row_visibilities if not x ]

        # Visible rows
        t_showing_rows_txt = 'Showing <sup id="{tid}_numrows" class="mqc_table_numrows">{nvisrows}</sup>/<sub>{nrows}</sub> rows'.format(tid=table_id, nvisrows=len(visible_rows), nrows=len(t_rows))

        # How many columns are visible?
        ncols_vis = (len(t_headers) + 1) - hidden_cols
        t_showing_cols_txt = ''
        if len(t_headers) > 1:
            t_showing_cols_txt = ' and <sup id="{tid}_numcols" class="mqc_table_numcols">{ncols_vis}</sup>/<sub>{ncols}</sub> columns'.format(tid=table_id, ncols_vis=ncols_vis, ncols=len(t_headers))

        # Build table header text
        html += """
        <small id="{tid}_numrows_text" class="mqc_table_numrows_text">{rows}{cols}.</small>
        """.format(tid=table_id, rows=t_showing_rows_txt, cols=t_showing_cols_txt)

    # Build the table itself
    collapse_class = 'mqc-table-collapse' if len(t_rows) > 10 and config.collapse_tables else ''
    html += """
        <div id="{tid}_container" class="mqc_table_container">
            <div class="table-responsive mqc-table-responsive {cc}">
                <table id="{tid}" class="table table-condensed mqc_table" data-title="{title}">
        """.format( tid=table_id, title=table_title, cc=collapse_class)

    # Build the header row
    col1_header = dt.pconfig.get('col1_header', 'Sample Name')
    html += '<thead><tr><th class="rowheader">{}</th>{}</tr></thead>'.format(col1_header, ''.join(t_headers.values()))

    # Build the table body
    html += '<tbody>'
    t_row_keys = t_rows.keys()
    if dt.pconfig.get('sortRows') is not False:
        t_row_keys = sorted(t_row_keys)
    for s_name in t_row_keys:
        # Hide the row if all cells are empty or hidden
        row_hidden = ' style="display:none"' if all(t_rows_empty[s_name].values()) else  ''
        html += '<tr{}>'.format(row_hidden)
        # Sample name row header
        html += '<th class="rowheader" data-original-sn="{sn}">{sn}</th>'.format(sn=s_name)
        for k in t_headers:
            html += t_rows[s_name].get(k, empty_cells[k])
        html += '</tr>'
    html += '</tbody></table></div>'
    if len(t_rows) > 10 and config.collapse_tables:
        html += '<div class="mqc-table-expand"><span class="glyphicon glyphicon-chevron-down" aria-hidden="true"></span></div>'
    html += '</div>'

    # Build the bootstrap modal to customise columns and order
    if not config.simple_output:
        html += """
    <!-- MultiQC Table Columns Modal -->
    <div class="modal fade" id="{tid}_configModal" tabindex="-1">
      <div class="modal-dialog modal-lg">
        <div class="modal-content">
          <div class="modal-header">
            <button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
            <h4 class="modal-title">{title}: Columns</h4>
          </div>
          <div class="modal-body">
            <p>Uncheck the tick box to hide columns. Click and drag the handle on the left to change order.</p>
            <p>
                <button class="btn btn-default btn-sm mqc_configModal_bulkVisible" data-target="#{tid}" data-action="showAll">Show All</button>
                <button class="btn btn-default btn-sm mqc_configModal_bulkVisible" data-target="#{tid}" data-action="showNone">Show None</button>
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
    </div> </div> </div>""".format( tid=table_id, title=table_title, trows=''.join(t_modal_headers.values()) )

    # Save the raw values to a file if requested
    if dt.pconfig.get('save_file') is True:
        fn = dt.pconfig.get('raw_data_fn', 'multiqc_{}'.format(table_id) )
        util_functions.write_data_file(dt.raw_vals, fn )
        report.saved_raw_data[fn] = dt.raw_vals

    return html
