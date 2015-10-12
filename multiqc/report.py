#!/usr/bin/env python

""" MultiQC report module. Holds the output from each
module. Is available to subsequent modules. """

from collections import defaultdict, OrderedDict

general_stats = {
    'headers': OrderedDict(),
    'rows': defaultdict(lambda:dict()),
    'shared_keys': defaultdict(lambda:dict())
}

