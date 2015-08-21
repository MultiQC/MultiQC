#!/usr/bin/env python

""" MultiQC modules base class, contains helper functions """

import json

class BaseMultiqcModule(object):

    def __init__(self):
        pass


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
                h = d[sn].keys()
                l.append(delim.join([''] + h))
            thesefields = [sn] + [ str(d[sn].get(k, '')) for k in h ]
            l.append( delim.join( thesefields ) )
        return '\n'.join(l)
