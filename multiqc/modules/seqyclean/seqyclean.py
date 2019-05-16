#!/usr/bin/env python3

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import linegraph, bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)
import traceback

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='My Module', anchor='mymod',
        href="http://www.awesome_bioinfo.com/my_module",
        info="is an example analysis module used for writing documentation.")

	 def parse_logs(self, f):
        data = {}
        for l in f.splitlines():
            s = l.split()
            data[s[0]] = s[1]
        return data
		
#print('Hello World')
