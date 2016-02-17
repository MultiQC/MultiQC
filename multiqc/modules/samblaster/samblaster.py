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
    """ Samblaster """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Samblaster', anchor='samblaster',
                                            href="https://github.com/GregoryFaust/samblaster",
                                            info="is a tools to remove duplicate molecules in a stream")

        self.samblaster_data = dict()
        for f in self.find_log_files(config.sp['samblaster'], filehandles=True):
            self.parse_samblaster(f)

        headers = OrderedDict()
        headers['n_tot'] = {
             'title': 'Read Pairs',
             'description': 'Total number of read pairs processed',
             'modify': lambda x: x / 1000000,
             'min': 0,
             'scale': 'RdYlGn',
             'format': '{:.0f}'
        }
        headers['n_dups'] = {
            'title': 'N Dups',
            'description': 'Number of duplicate reads filtered by samblaster',
            'modify': lambda x: x / 1000000,
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:.0f}'
        }
        headers['pct_dups'] = {
            'title': 'Pct Dups',
            'description': 'Percent Duplicates',
            'scale': 'RdYlGn-rev',
            'max': 100,
            'min': 0
        }

        self.general_stats_addcols(self.samblaster_data, headers)

        # Write parsed report data to a file
        self.write_data_file(self.samblaster_data, 'multiqc_samblaster')

        if len(self.samblaster_data) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.samblaster_data)))

    def parse_samblaster(self, f):
        """ Go through log file looking for samblaster output """
        dups_regex = "samblaster: Removed (\d+) of (\d+) \((\d+.\d+)%\) read ids as duplicates"
        name_regex = "@RG\\\\tID:(\S*?)\\\\t"
        data = {}
        fh = f['f']
        for l in fh:
            match = re.search(name_regex, l)
            if match:
                data['s_name'] = match.group(1)

            match = re.search(dups_regex, l)
            if match:
                data['n_dups'] = match.group(1)
                data['n_tot'] = match.group(2)
                data['pct_dups'] = match.group(3)

        if 's_name' in data:
            s_name = data['s_name']
            self.add_data_source(f, s_name)
            self.samblaster_data[s_name] = dict(n_dups=int(data['n_dups']),
                                                n_tot=int(data['n_tot']),
                                                pct_dups=float(data['pct_dups']))

    def add_skewer_data(self, s_name, data, f):
        pass
        # stats = ['r_processed', 'r_short_filtered', 'r_empty_filtered', 'r_avail', 'r_trimmed', 'r_untrimmed']
        # if s_name in self.skewer_data:
        #     log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
        # self.skewer_data[s_name] = {}
        # self.add_data_source(f, s_name)
        # for k in stats:
        #     self.skewer_data[s_name][k] = int(data[k])
        #
        # self.skewer_data[s_name]['pct_avail'] = 100.0 * float(data['r_avail']) / float(data['r_processed'])
        # self.skewer_data[s_name]['pct_trimmed'] = 100.0 * float(data['r_trimmed']) / float(data['r_avail'])
        # self.skewer_data[s_name]['pct_untrimmed'] = 100.0 * float(data['r_untrimmed']) / float(data['r_avail'])
