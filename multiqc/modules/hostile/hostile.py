#!/usr/bin/env python

import json
import logging
import os
from collections import OrderedDict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """Hostile Module"""
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='Hostile',
            anchor='hostile',
            href='https://github.com/bede/hostile',
            info='Short and Long host reads removal tools.',
            doi="10.1093/bioinformatics/btad728",
        )

        self.data = dict()
        for f in self.find_log_files('hostile', filehandles=True):
            self.parse_logs(f)
            

        if len(self.data) == 0:
            log.debug('No reports found: hostile')
            raise UserWarning

        self.hostile_plot()

    def parse_logs(self, json_file):
        """
            Parsing json_file
        """
        try:
            data = json.load(json_file['f'])
                
        except json.JSONDecodeError as e:
            log.warning(f"Could not parse JSON file {json_file['f']}")
            return

        s_name = self.clean_s_name(json_file['fn'])
        self.data[s_name] = data
        #print(self.data)

    def hostile_plot(self):
        """
        Extract the data
        """
        data = {}
        for f_name, values in self.data.items():
            s_name = values[0]['fastq1_in_name'].split(".")[0]
            database = os.path.basename(values[0]['index'])
            data[s_name] = {
                'Cleaned reads': values[0]['reads_out'],
                'Host reads': values[0]['reads_removed']
            }
        cat = ['Cleaned reads', 'Host reads']

        pconfig = {
            'title': 'Hostile: Reads Filtered',
            'plot_type': 'bargraph',
            'stacked_data': True
        }

        self.add_section(
            name='Reads Filtering',
            anchor='hostile-reads',
            description=f'This plot shows the number of cleaned reads vs host-reads per sample (database index: {database}).',
            plot=bargraph.plot(data, cat, pconfig)
        )