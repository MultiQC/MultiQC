#!/usr/bin/env python

""" MultiQC Utility functions, used in a variety of places. """

from __future__ import print_function
import io
import json
import os
import yaml

from multiqc import config

def write_data_file(data, fn, sort_cols=False, data_format=None):
    """ Write a data file to the report directory. Will not do anything
    if config.data_dir is not set.
    :param: data - a 2D dict, first key sample name (row header),
            second key field (column header). 
    :param: fn - Desired filename. Directory will be prepended automatically.
    :param: sort_cols - Sort columns alphabetically
    :param: data_format - Output format. Defaults to config.data_format (usually tsv)
    :return: None """
    
    if config.data_dir is not None:
        
        # Add relevant file extension to filename
        if data_format is None:
            data_format = config.data_format
        fn = '{}.{}'.format(fn, config.data_format_extensions[data_format])
        
        # Save file
        with io.open (os.path.join(config.data_dir, fn), 'w', encoding='utf-8') as f:
            if data_format == 'json':
                jsonstr = json.dumps(data, indent=4, ensure_ascii=False)
                print( jsonstr.encode('utf-8', 'ignore').decode('utf-8'), file=f)
            elif data_format == 'yaml':
                yaml.dump(data, f, default_flow_style=False)
            else:
                # Default - tab separated output
                # Get all headers
                h = ['Sample']
                for sn in sorted(data.keys()):
                    for k in data[sn].keys():
                        if type(data[sn][k]) is not dict and k not in h:
                            h.append(str(k))
                if sort_cols:
                    h = sorted(h)
                
                # Get the rows
                rows = [ "\t".join(h) ]
                for sn in sorted(data.keys()):
                    # Make a list starting with the sample name, then each field in order of the header cols
                    l = [sn] + [ str(data[sn].get(k, '')) for k in h[1:] ]
                    rows.append( "\t".join(l) )
                
                body = '\n'.join(rows)
                
                print( body.encode('utf-8', 'ignore').decode('utf-8'), file=f)
