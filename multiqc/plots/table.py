#!/usr/bin/env python

""" MultiQC functions to plot a report beeswarm group """

import json
import logging
import os
import random

from multiqc.utils import report, config
logger = logging.getLogger(__name__)

letters = 'abcdefghijklmnopqrstuvwxyz'

def plot (data, headers=[], pconfig={}):
    """ Return HTML for a MultiQC table. See CONTRIBUTING.md for
    further instructions on use.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param headers: optional dict with column config in key:value pairs.
    :return: HTML ready to be inserted into the page
    """
    
    return make_table_html( datatable(data, headers, pconfig) )
    
    

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
        for idx, d in data.items():
            
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
                    for (sname, samp) in data.items():
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
        for idx in headers.keys():
            for k in headers.keys():
                sk = headers[idx][k].get('shared_key', None)
                if sk is not None:
                    shared_keys[sk]['scale'] = headers[idx][k]['scale']
                    shared_keys[sk]['dmax']  = max(headers[idx][k]['dmax'], shared_keys[sk].get('dmax', headers[idx][k]['dmax']))
                    shared_keys[sk]['dmin']  = max(headers[idx][k]['dmin'], shared_keys[sk].get('dmin', headers[idx][k]['dmin']))
        
        # Overwrite shared key settings
        for idx in headers.keys():
            for k in headers[idx].keys():
                sk = headers[idx][k].get('shared_key', None)
                if sk is not None:
                    headers[idx][k]['scale'] = shared_keys[sk]['scale']
                    headers[idx][k]['dmax'] = shared_keys[sk]['dmax']
                    headers[idx][k]['dmin'] = shared_keys[sk]['dmin']
        
        # Assign to class
        self.data = data
        self.headers = headers
        self.pconfig = pconfig


def make_table_html (dt):
    """
    Build the HTML needed for a MultiQC table.
    :param data: MultiQC datatable object
    """
    
    for idx, d in dt.data.items():
        for k in d.keys():
            
            rid = headers[idx][k]['rid']
            
            if headers[idx][k].get('shared_key', None) is not None:
                sk = ' data-shared-key={}'.format(sk)
            else:
                sk = ''
            
            general_stats_html['headers'][rid] = '<th \
                    id="header_{rid}" \
                    class="chroma-col {rid}" \
                    data-chroma-scale="{scale}" \
                    data-chroma-max="{max}" \
                    data-chroma-min="{min}" \
                    {sk}><span data-toggle="tooltip" title="{mod}: {descrip}">{title}</span></th>' \
                        .format(rid=rid, scale=headers[idx][k]['scale'], max=headers[idx][k]['dmax'],
                            min=headers[idx][k]['dmin'], sk=sk, mod=mod,
                            descrip=headers[idx][k]['description'], title=headers[idx][k]['title'])
            
            # Add the data table cells
            nrows = 0
            for (sname, samp) in dt.data[idx]['data'].items():
                if k in samp:
                    val = samp[k]
                    general_stats_raw[sname][rid] = val
                    
                    if 'modify' in headers[idx][k] and callable(headers[idx][k]['modify']):
                        val = headers[idx][k]['modify'](val)
                    
                    try:
                        dmin = headers[idx][k]['dmin']
                        dmax = headers[idx][k]['dmax']
                        percentage = ((float(val) - dmin) / (dmax - dmin)) * 100;
                        percentage = min(percentage, 100)
                        percentage = max(percentage, 0)
                    except (ZeroDivisionError,ValueError):
                        percentage = 0
                    
                    try:
                        val = headers[idx][k]['format'].format(val)
                    except ValueError:
                        try:
                            val = headers[idx][k]['format'].format(float(samp[k]))
                        except ValueError:
                            val = samp[k]
                    except:
                        val = samp[k]
                    
                    general_stats_html['rows'][sname][rid] = \
                        '<td class="data-coloured {rid}" >\
                            <div class="wrapper">\
                                <span class="bar" style="width:{percentage}%;"></span>\
                                <span class="val">{val}</span>\
                            </div>\
                        </td>'.format(rid=rid, percentage=percentage, val=val)
                    nrows += 1
            
            # Remove header if we don't have any filled cells for it
            if nrows == 0:
                general_stats_html['headers'].pop(rid, None)
                logger.debug('Removing header {} from general stats table, as no data'.format(k))
    
    
    
    
    
    