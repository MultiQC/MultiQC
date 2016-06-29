#!/usr/bin/env python

""" MultiQC functions to plot a table """

from collections import defaultdict, OrderedDict
import json
import logging
import os
import random

from multiqc.utils import config, util_functions
from multiqc.plots import table_object, beeswarm
logger = logging.getLogger(__name__)

letters = 'abcdefghijklmnopqrstuvwxyz'

def plot (data, headers=[], pconfig={}):
    """ Return HTML for a MultiQC table.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param headers: list of optional dicts with column config in key:value pairs.
    :return: HTML ready to be inserted into the page
    """
    
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
        warning = """<p>This report is very large.
            A beeswarm plot has been generated instead of a table to aid usability.
            To disable this, set <code>max_table_rows</code> to a very high number
            in your MultiQC config file.</p>"""
        return warning + beeswarm.make_plot( dt )
    else:
        return make_table ( dt )


def make_table (dt):
    """
    Build the HTML needed for a MultiQC table.
    :param data: MultiQC datatable object
    """
    
    table_id = dt.pconfig.get('id', 'table_{}'.format(''.join(random.sample(letters, 4))) )
    t_headers = OrderedDict()
    t_modal_headers = OrderedDict()
    t_rows = defaultdict(lambda: dict())
    dt.raw_vals = defaultdict(lambda: dict())
    empty_cells = dict()
    hidden_cols = 1
    
    for idx, hs in enumerate(dt.headers):
        for k, header in hs.items():
            
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
            
            data_attr = 'data-chroma-scale="{}" data-chroma-max="{}" data-chroma-min="{}" {}' \
                .format(header['scale'], header['dmax'], header['dmin'], shared_key)
            
            cell_contents = '<span data-toggle="tooltip" title="{}: {}">{}</span>' \
                .format(header['namespace'], header['description'], header['title'])
            
            t_headers[rid] = '<th id="header_{rid}" class="chroma-col {rid} {h}" {d}>{c}</th>' \
                .format(rid=rid, d=data_attr, h=hide, c=cell_contents)
            
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
            
            # Add the data table cells
            for (s_name, samp) in dt.data[idx].items():
                if k in samp:
                    val = samp[k]
                    dt.raw_vals[s_name][rid] = val
                    
                    if 'modify' in header and callable(header['modify']):
                        val = header['modify'](val)
                    
                    try:
                        dmin = header['dmin']
                        dmax = header['dmax']
                        percentage = ((float(val) - dmin) / (dmax - dmin)) * 100;
                        percentage = min(percentage, 100)
                        percentage = max(percentage, 0)
                    except (ZeroDivisionError,ValueError):
                        percentage = 0
                    
                    try:
                        val = header['format'].format(val)
                    except ValueError:
                        try:
                            val = header['format'].format(float(samp[k]))
                        except ValueError:
                            val = samp[k]
                    except:
                        val = samp[k]
                    
                    # Build HTML
                    bar_html = '<span class="bar" style="width:{}%;"></span>'.format(percentage)
                    val_html = '<span class="val">{}</span>'.format(val)
                    wrapper_html = '<div class="wrapper">{}{}</div>'.format(bar_html, val_html)
                    
                    t_rows[s_name][rid] = \
                        '<td class="data-coloured {rid} {h}">{c}</td>'.format(rid=rid, h=hide, c=wrapper_html)
            
            # Remove header if we don't have any filled cells for it
            if sum([len(rows) for rows in t_rows.values()]) == 0:
                t_headers.pop(rid, None)
                t_modal_headers.pop(rid, None)
                logger.debug('Removing header {} from general stats table, as no data'.format(k))
    
    #
    # Put everything together
    #
    
    # Buttons above the table
    html = """
        <button type="button" class="mqc_table_copy_btn btn btn-default btn-sm" data-clipboard-target="#{tid}">
            <span class="glyphicon glyphicon-copy"></span> Copy table
        </button>
        <button type="button" class="mqc_table_configModal_btn btn btn-default btn-sm" data-toggle="modal" data-target="#{tid}_configModal">
            <span class="glyphicon glyphicon-th"></span> Configure Columns
        </button>
        <button type="button" class="mqc_table_sortHighlight btn btn-default btn-sm" data-target="#{tid}" data-direction="desc" style="display:none;">
            <span class="glyphicon glyphicon-sort-by-attributes-alt"></span> Sort by highlight
        </button>
        <small id="{tid}_numrows_text" class="mqc_table_numrows_text">Showing <sup id="{tid}_numrows" class="mqc_table_numrows">{nrows}</sup>/<sub>{nrows}</sub> rows and <sup id="{tid}_numcols" class="mqc_table_numcols">{ncols_vis}</sup>/<sub>{ncols}</sub> columns.</small>
    """.format(tid=table_id, nrows=len(t_rows), ncols_vis = len(t_headers)-hidden_cols, ncols=len(t_headers))
    
    # Build the table itself
    html += """
        <div id="{tid}_container" class="mqc_table_container">
            <div class="table-responsive">
                <table id="{tid}" class="table table-condensed mqc_table">
        """.format(tid=table_id)
    
    # Build the header row
    html += '<thead><tr><th class="rowheader">Sample Name</th>{}</tr></thead>'.format(''.join(t_headers.values()))
    
    # Build the table body
    html += '<tbody>'
    for s_name in sorted(t_rows.keys()):
        html += '<tr>'
        # Sample name row header
        html += '<th class="rowheader" data-original-sn="{sn}">{sn}</th>'.format(sn=s_name)
        for k in t_headers:
            html += t_rows[s_name].get(k, empty_cells[k])
        html += '</tr>'
    html += '</tbody></table></div></div>'
    
    # Build the bootstrap modal to customise columns and order
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
            <table class="table mqc_table mqc_sortable mqc_configModal_table" id="{tid}_configModal_table">
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
    </div> </div> </div>""".format( tid=table_id, title=dt.pconfig.get('table_title', table_id), trows=''.join(t_modal_headers.values()) )
    
    # Save the raw values to a file if requested
    if dt.pconfig.get('save_file') is True:
        fn = dt.pconfig.get('raw_data_fn', 'multiqc_{}'.format(table_id) )
        util_functions.write_data_file(dt.raw_vals, fn )
    
    return html
    
    
    
    
    
    