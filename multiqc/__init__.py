#!/usr/bin/env python

""" MultiQC modules base class, contains helper functions """

import json
import os

class BaseMultiqcModule(object):

    def __init__(self):
        pass

    def clean_s_name(self, s_name, root, prepend_dirs=False, trimmed=True):
        """ Helper function to take a long file name and strip it
        back to a clean sample name. Somewhat arbitrary. """
        # Split then take first section to remove everything after these matches
        s_name = s_name.split(".gz",1)[0]
        s_name = s_name.split(".fastq",1)[0]
        s_name = s_name.split(".fq",1)[0]
        s_name = s_name.split(".bam",1)[0]
        s_name = s_name.split(".sam",1)[0]
        if trimmed:
            s_name = s_name.split("_val_1",1)[0]
            s_name = s_name.split("_val_2",1)[0]
            s_name = s_name.split("_trimmed",1)[0]
        if prepend_dirs:
            s_name = "{} | {}".format(root.replace(os.sep, ' | '), s_name).lstrip('. | ')
        return s_name
    
    def dict_to_csv (self, d, delim="\t"):
        """ Takes a 2D dictionary and returns a string suitable for
        writing to a .csv file. First key should be sample name
        (row header), second key should be field (column header)
        Takes dictionary as input, optional 'delim' field to specify
        column delimiter (default: tab) """

        h = None # We make a list of keys to ensure consistent order
        l = list()
        for sn in d:
            if h is None:
                h = list(d[sn].keys())
                l.append(delim.join([''] + h))
            thesefields = [sn] + [ str(d[sn].get(k, '')) for k in h ]
            l.append( delim.join( thesefields ) )
        return ('\n'.join(l)).encode('utf-8', 'ignore').decode('utf-8')
