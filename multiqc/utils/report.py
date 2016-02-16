#!/usr/bin/env python

""" MultiQC report module. Holds the output from each
module. Is available to subsequent modules. Contains
helper functions to generate markup for report. """

from collections import defaultdict, OrderedDict
import logging

import multiqc
from multiqc import (logger)
from multiqc.utils.log import init_log, LEVELS

general_stats = OrderedDict()
general_stats_html = {
    'headers': OrderedDict(),
    'rows': defaultdict(lambda:dict())
}
general_stats_raw = defaultdict(lambda:OrderedDict())
data_sources = defaultdict(lambda:defaultdict(lambda:defaultdict()))



def general_stats_build_table():
    """ Helper function to add to the General Statistics table.
    Parses report.general_stats and returns HTML for general stats table.
    Also creates report.general_stats_raw for multiqc_general_stats.txt
    :param data: A dict with the data. First key should be sample name,
                 then the data key, then the data.
    :param headers: Dict / OrderedDict with information for the headers, 
                    such as colour scales, min and max values etc.
                    See docs/writing_python.md for more information.
    :return: None
    """
    
    # First - collect settings for shared keys
    shared_keys = defaultdict(lambda: dict())
    for mod in general_stats.keys():
        headers = general_stats[mod]['headers']
        for k in headers.keys():
            sk = headers[k].get('shared_key', None)
            if sk is not None:
                shared_keys[sk]['scale'] = headers[k]['scale']
                shared_keys[sk]['dmax']  = max(headers[k]['dmax'], shared_keys[sk].get('dmax', headers[k]['dmax']))
                shared_keys[sk]['dmin']  = max(headers[k]['dmin'], shared_keys[sk].get('dmin', headers[k]['dmin']))
    
    # Now build required HTML
    modcols = ['228,26,28', '55,126,184', '77,175,74', '152,78,163', '255,127,0', '255,255,51', '166,86,40', '247,129,191', '153,153,153']
    midx = 0
    for mod in general_stats.keys():
        headers = general_stats[mod]['headers']
        for k in headers.keys():
            
            rid = headers[k]['rid']
            
            headers[k]['modcol'] = modcols[midx]
            
            # Overwrite config with shared key settings
            sk = headers[k].get('shared_key', None)
            if sk is not None:
                headers[k]['scale'] = shared_keys[sk]['scale']
                headers[k]['dmax'] = shared_keys[sk]['dmax']
                headers[k]['dmin'] = shared_keys[sk]['dmin']
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
                        .format(rid=rid, scale=headers[k]['scale'], max=headers[k]['dmax'],
                            min=headers[k]['dmin'], sk=sk, mod=mod,
                            descrip=headers[k]['description'], title=headers[k]['title'])
            
            # Add the data table cells
            nrows = 0
            for (sname, samp) in general_stats[mod]['data'].items():
                if k in samp:
                    val = samp[k]
                    general_stats_raw[sname][rid] = val
                    
                    if 'modify' in headers[k] and callable(headers[k]['modify']):
                        val = headers[k]['modify'](val)
                    
                    try:
                        dmin = headers[k]['dmin']
                        dmax = headers[k]['dmax']
                        percentage = ((float(val) - dmin) / (dmax - dmin)) * 100;
                        percentage = min(percentage, 100)
                        percentage = max(percentage, 0)
                    except (ZeroDivisionError,ValueError):
                        percentage = 0
                    
                    try:
                        val = headers[k]['format'].format(val)
                    except ValueError:
                        try:
                            val = headers[k]['format'].format(float(samp[k]))
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
        
        # Index for colouring by module
        midx += 1
        if midx > (len(modcols) - 1):
            midx = 0
    
    return None
    
    
    
def dict_to_csv (d, delim="\t", sort_cols=False):
    """ Converts a dict to a CSV string
    :param d: 2D dictionary, first keys sample names and second key
              column headers. If second key is not a string, it is skipped.
    :param delim: optional delimiter character. Default: \t
    :return: Flattened string, suitable to write to a CSV file.
    """
    
    # Get all headers
    h = ['Sample']
    for sn in sorted(d.keys()):
        for k in d[sn].keys():
            if type(d[sn][k]) is not dict and k not in h:
                h.append(k)
    if sort_cols:
        h = sorted(h)
    
    # Get the rows
    rows = [ delim.join(h) ]
    for sn in sorted(d.keys()):
        # Make a list starting with the sample name, then each field in order of the header cols
        l = [sn] + [ str(d[sn].get(k, '')) for k in h[1:] ]
        rows.append( delim.join(l) )
    
    body = '\n'.join(rows)
    
    return body.encode('utf-8', 'ignore').decode('utf-8')

def data_sources_to_csv (d, delim="\t"):
    lines = [['Module', 'Section', 'Sample Name', 'Source']]
    for mod in d:
        for sec in d[mod]:
            for s_name, source in d[mod][sec].items():
                lines.append([mod, sec, s_name, source])
    body = '\n'.join([delim.join(l) for l in lines])
    return body.encode('utf-8', 'ignore').decode('utf-8')
    
    