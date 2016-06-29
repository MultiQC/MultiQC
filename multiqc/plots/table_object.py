#!/usr/bin/env python

""" MultiQC datatable class, used by tables and beeswarm plots """

from collections import defaultdict
import logging
import random

from multiqc.utils import config

logger = logging.getLogger(__name__)

letters = 'abcdefghijklmnopqrstuvwxyz'

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
            except (IndexError, AssertionError):
                keys = d.keys()
            
            # Check that we have some data in each column
            empties = list()
            for k in keys:
                n = 0
                for samp in data[idx].values():
                    if k in samp:
                        n += 1
                if n == 0:
                    empties.append(k)
            for k in empties:
                keys = [j for j in keys if j != k]
                del headers[idx][k]
            
            for k in keys:
                # Unique id to avoid overwriting by other datasets
                headers[idx][k]['rid'] = '{}_{}'.format(''.join(random.sample(letters, 4)), k)
                
                # Use defaults / data keys if headers not given
                headers[idx][k]['namespace']   = headers[idx][k].get('namespace', '')
                headers[idx][k]['title']       = headers[idx][k].get('title', k)
                headers[idx][k]['description'] = headers[idx][k].get('description', headers[idx][k]['title'])
                headers[idx][k]['scale']       = headers[idx][k].get('scale', 'GnBu')
                headers[idx][k]['format']      = headers[idx][k].get('format', '{:.1f}')
                if 'colour' not in headers[idx][k]:
                    cidx = idx
                    while cidx >= len(sectcols):
                        cidx -= len(sectcols)
                    headers[idx][k]['colour'] = sectcols[cidx]
                
                # Overwrite hidden if set in user config
                try:
                    # Config has True = visibile, False = Hidden. Here we're setting "hidden" which is inverse
                    headers[idx][k]['hidden'] = not config.table_columns_visible[ headers[idx][k]['namespace'] ][k]
                except KeyError:
                    pass
                
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
                    for s_name, samp in data[idx].items():
                        try:
                            val = float(samp[k])
                            if 'modify' in headers[idx][k] and callable(headers[idx][k]['modify']):
                                val = float(headers[idx][k]['modify'](val))
                            if setdmax:
                                headers[idx][k]['dmax'] = max(headers[idx][k]['dmax'], val)
                            if setdmin:
                                headers[idx][k]['dmin'] = min(headers[idx][k]['dmin'], val)
                        except ValueError:
                            val = samp[k] # couldn't convert to float - keep as a string
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
