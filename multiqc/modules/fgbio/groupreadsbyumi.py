#!/usr/bin/env python

""" MultiQC module to parse output from fgbio GroupReadsByUmi
"""

############################################################
######  LOOKING FOR AN EXAMPLE OF HOW MULTIQC WORKS?  ######
############################################################
#### Stop! This is one of the most complicated modules. ####
#### Have a look at Kallisto for a simpler example.     ####
############################################################

from __future__ import print_function
from collections import OrderedDict
import io
import json
import logging
import os
import re
import zipfile

from multiqc import config
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Fgbio', anchor='fgbio',
        href="http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html",
        info="fgbio is a command line toolkit for working with genomic and particularly next generation sequencing data.")

        self.umi_data = dict()

        # Find and parse unzipped FastQC reports
        for f in self.find_log_files('fgbio/groupreadsbyumi'):
            self.parse_hist_report(f)

        log.info("Processed {} report(s)".format(len(self.umi_data)))

        self.add_section(name = 'GroupReadsByUmi statistics',
        anchor = 'fgbio-groupreadsbyumi',
        description = 'During GroupReadsByUmi processing family size count data is generated, showing number of UMIs represented by a certain number of reads. ',
        helptext = '''**Note!** Multiqc expects the input file to have the following format <SAMPLE>.histo.tsv, SAMPLE will be for naming data.
                    ''',
        plot = self.create_plot())


    def parse_hist_report(self, f):
        sample_name = self.get_s_name(f)
        family_size = []
        for line in f['f'].splitlines():
            if not line.startswith("family_size"):
                family_size.append(tuple(line.split("\t")))

        self.umi_data[sample_name] = { int(s):int(d[1]) for s, d in enumerate(family_size,1)}


    def create_plot(self):
        config = {
            'id': 'umi_support',
            'title': 'Familizy size count',
            'ylab': '# UMI',
            'xlab': 'Reads supporting UMI',
            "xmax": 15
        }

        return linegraph.plot(self.umi_data, config)

    def get_s_name(self, f):
        sample_name = f['s_name']
        if sample_name.endswith('.histo'):
            sample_name = sample_name[:-6]
        return sample_name
