#!/usr/bin/env python

""" MultiQC functions to plot a table """

from collections import defaultdict
import json
import logging
import os
import random

from multiqc.utils import config
logger = logging.getLogger(__name__)

letters = 'abcdefghijklmnopqrstuvwxyz'

def plot (data, headers=[], pconfig={}):
    """ Return HTML for a MultiQC table.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param headers: list of optional dicts with column config in key:value pairs.
    :return: HTML ready to be inserted into the page
    """
    
    return make_table ( datatable(data, headers, pconfig) )
    
    

class datatable (object):
    
    def __init__ (self, data, headers=[], pconfig={}):
        """ Prepare data for use in a table or plot """
        
        # Given one dataset - turn it into a list
        if type(data) is not list:
            data = [data]
        if type(headers) is not list:
            headers = [headers]
        
        sectcols = ['55,126,184', '77,175,74', '152,78,163', '255,127,0', '228,26,28', '255,255,51', '166,86,40', '247,129,191', '153,153,153']
        shared_keys = defaultdict(lambda: dict())
        
        # Go through each table section
        for idx, d in enumerate(data):
            
            # Get the header keys
            try:
                keys = headers[idx].keys()
                assert len(keys) > 0
            except IndexError, AssertionError:
                keys = d.keys()
            
            for k in keys:
                # Unique id to avoid overwriting by other datasets
                headers[idx][k]['rid'] = '{}_{}'.format(''.join(random.sample(letters, 4)), k)
                
                # Use defaults / data keys if headers not given
                headers[idx][k]['title']       = headers[idx][k].get('title', k)
                headers[idx][k]['description'] = headers[idx][k].get('description', headers[idx][k]['title'])
                headers[idx][k]['scale']       = headers[idx][k].get('scale', 'GnBu')
                headers[idx][k]['format']      = headers[idx][k].get('format', '{:.1f}')
                if 'colour' not in headers[idx][k]:
                    cidx = idx
                    while cidx >= len(sectcols):
                        cidx -= len(sectcols)
                    headers[idx][k]['colour'] = sectcols[cidx]
                
                # Work out max and min value if not given
                setdmax = False
                setdmin = False
                try:
                    headers[idx][k]['dmax'] = float(headers[idx][k]['max'])
                except KeyError:
                    headers[idx][k]['dmax'] = float("-inf")
                    setdmax = True
                
                try:
                    headers[idx][k]['dmin'] = float(headers[idx][k]['min'])
                except KeyError:
                    headers[idx][k]['dmin'] = float("inf")
                    setdmin = True
                
                # Figure out the min / max if not supplied
                if setdmax or setdmin:
                    for (s_name, samp) in enumerate(data):
                        try:
                            val = float(samp[k])
                            if 'modify' in headers[idx][k] and callable(headers[idx][k]['modify']):
                                val = float(headers[idx][k]['modify'](val))
                            if setdmax:
                                headers[idx][k]['dmax'] = max(headers[idx][k]['dmax'], val)
                            if setdmin:
                                headers[idx][k]['dmin'] = min(headers[idx][k]['dmin'], val)
                        except KeyError:
                            pass # missing data - skip
        
        # Collect settings for shared keys
        shared_keys = defaultdict(lambda: dict())
        for idx, hs in enumerate(headers):
            for k in hs.keys():
                sk = headers[idx][k].get('shared_key', None)
                if sk is not None:
                    shared_keys[sk]['scale'] = headers[idx][k]['scale']
                    shared_keys[sk]['dmax']  = max(headers[idx][k]['dmax'], shared_keys[sk].get('dmax', headers[idx][k]['dmax']))
                    shared_keys[sk]['dmin']  = max(headers[idx][k]['dmin'], shared_keys[sk].get('dmin', headers[idx][k]['dmin']))
        
        # Overwrite shared key settings
        for idx, hs in enumerate(headers):
            for k in hs.keys():
                sk = headers[idx][k].get('shared_key', None)
                if sk is not None:
                    headers[idx][k]['scale'] = shared_keys[sk]['scale']
                    headers[idx][k]['dmax'] = shared_keys[sk]['dmax']
                    headers[idx][k]['dmin'] = shared_keys[sk]['dmin']
        
        # Assign to class
        self.data = data
        self.headers = headers
        self.pconfig = pconfig


def make_table (dt):
    """
    Build the HTML needed for a MultiQC table.
    :param data: MultiQC datatable object
    """
    
    t_headers = dict()
    t_rows = defaultdict(lambda: dict())
    raw_vals = defaultdict(lambda: dict())
    
    for idx, hs in enumerate(dt.headers):
        for k, header in hs.items():
            
            rid = header['rid']
            
            # Build the table header cell
            if header.get('shared_key', None) is not None:
                shared_key = ' data-shared-key={}'.format(header['shared_key'])
            else:
                shared_key = ''
            
            data_attr = 'data-chroma-scale="{}" data-chroma-max="{}" data-chroma-min="{}" {}' \
                .format(header['scale'], header['dmax'], header['dmin'], shared_key)
            
            cell_contents = '<span data-toggle="tooltip" title="{}">{}</span>' \
                .format(header['description'], header['title'])
            
            t_headers[rid] = '<th id="header_{rid}" class="chroma-col {rid}" {d}>{c}</th>' \
                .format(rid=rid, d=data_attr, c=cell_contents)
            
            # Add the data table cells
            for (s_name, samp) in dt.data[idx].items():
                if k in samp:
                    val = samp[k]
                    raw_vals[s_name][rid] = val
                    
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
                        '<td class="data-coloured {rid}" >{c}</td>'.format(rid=rid, c=wrapper_html)
            
            # Remove header if we don't have any filled cells for it
            if sum([len(rows) for rows in t_rows.values()]) == 0:
                t_headers.pop(rid, None)
                logger.debug('Removing header {} from general stats table, as no data'.format(k))
    
    
    # Put everything together into a table
    html = '<div class="table-responsive">'
    html += '<table id="{}" class="table table-condensed mqc_table">'.format(dt.pconfig.get('id', ''))
    
    # Build the header row
    html += '<thead><tr><th class="rowheader">Sample Name</th>{}</tr></thead>'.format(''.join(t_headers))
    
    # Build the table body
    html += '<tbody>'
    for s_name in sorted(t_rows.keys()):
        html += '<tr>'
        # Sample name row header
        html += '<th class="rowheader" data-original-sn="{sn}">{sn}</th>'.format(sn=s_name)
        for k in t_headers:
            html += t_rows[s_name].get(k, '<td class="{}"></td>'.format(k) )
        html += '</tr>'
    html += '</tbody></table></div>'
    
    return html
    
    
    
    
    
    