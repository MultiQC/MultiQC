#!/usr/bin/env python

""" MultiQC module to parse output from fgbio GroupReadsByUmi
"""

from __future__ import print_function
from collections import OrderedDict

import logging

from multiqc import config
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class GroupReadsByUmiMixin():

    def parse_groupreadsbyumi_log(self):
        umi_data = dict()

        for f in self.find_log_files('fgbio/groupreadsbyumi'):
            sample_name = f['s_name']
            family_size = []
            for line in f['f'].splitlines():
                if not line.startswith("family_size"):
                    family_size.append(tuple(line.split("\t")))

            umi_data[sample_name] = { int(s):int(d[1]) for s, d in enumerate(family_size,1)}

        if not umi_data:
            raise UserWarning("No fgbio GroupReadsByUmi logs found")

        log.info("Processed {} report(s)".format(len(self.umi_data)))

        return umi_data
