#!/usr/bin/env python

""" MultiQC module to parse output from QualiMap """

from __future__ import print_function

import os
from collections import OrderedDict
import logging
import re
from multiqc import config, BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """ Skewer """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Skewer', anchor='skewer',
                                            href="https://github.com/relipmoc/skewer.git",
                                            info="is a tools to trim adapters off of reads")

        self.skewer_data = dict()
        for f in self.find_log_files(config.sp['skewer'], filehandles=True):
            self.parse_skewer_log(f)

        headers = OrderedDict()
        headers['r_processed'] = {
            'title': 'Read Pairs',
            'description': 'Total number of read pairs processed',
            'modify': lambda x: x / 1000000,
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:.0f}'
        }
        headers['r_avail'] = {
            'title': 'Kept',
            'description': 'Read pairs available after trimming',
            'modify': lambda x: x / 1000000,
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:.0f}'
        }
        headers['pct_trimmed'] = {
            'title': 'Pct Trimmed',
            'description': 'Percent Trimmed',
            'scale': 'RdYlGn-rev',
            'max': 100,
            'min': 0
        }

        self.general_stats_addcols(self.skewer_data, headers)

        # Write parsed report data to a file
        self.write_data_file(self.skewer_data, 'multiqc_skewer')

        if len(self.skewer_data) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning

        log.debug("Found {} reports".format(len(self.skewer_data)))

    def parse_skewer_log(self, f):
        """ Go through log file looking for skewer output """
        fh = f['f']
        regexes = {
            'fq1': "Input file:\s+(\S+).(fastq$|fastq.gz$)",
            'fq2': "Paired file:\s+(\S+).(fastq$|fastq.gz$)",
            'r_processed': "(\d+) read|reads pairs? processed",
            'r_short_filtered': "(\d+) \(\s*\d+.\d+%\) short read",
            'r_empty_filtered': "(\d+) \(\s*\d+.\d+%\) empty read",
            'r_avail': "(\d+) \(\s*\d+.\d+%\) read",
            'r_trimmed': "(\d+) \(\s*\d+.\d+%\) trimmed read",
            'r_untrimmed': "(\d+) \(\s*\d+.\d+%\) untrimmed read"
        }

        data = dict()
        for k, v in regexes.items():
            data[k] = 0
        data['fq1'] = None
        data['fq2'] = None

        for l in fh:
            for k, r in regexes.items():
                match = re.search(r, l)
                if match:
                    data[k] = match.group(1).replace(',', '')

        if data['fq1'] is not None:
            s_name = os.path.basename(data['fq1'])
            self.add_skewer_data(s_name, data, f)

        if data['fq2'] is not None:
            s_name = os.path.basename(data['fq2'])
            self.add_skewer_data(s_name, data, f)

    def add_skewer_data(self, s_name, data, f):
        stats = ['r_processed', 'r_short_filtered', 'r_empty_filtered', 'r_avail', 'r_trimmed', 'r_untrimmed']
        if s_name in self.skewer_data:
            log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
        self.skewer_data[s_name] = {}
        self.add_data_source(f, s_name)
        for k in stats:
            self.skewer_data[s_name][k] = int(data[k])

        self.skewer_data[s_name]['pct_avail'] = 100.0 * float(data['r_avail']) / float(data['r_processed'])
        self.skewer_data[s_name]['pct_trimmed'] = 100.0 * float(data['r_trimmed']) / float(data['r_avail'])
        self.skewer_data[s_name]['pct_untrimmed'] = 100.0 * float(data['r_untrimmed']) / float(data['r_avail'])

