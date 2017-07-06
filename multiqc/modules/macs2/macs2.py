#!/usr/bin/env python

""" MultiQC module to parse output from MACS2 """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='MACS2', anchor='macs',
        href='https://github.com/taoliu/MACS',
        info="identifies transcription factor binding sites in ChIP-seq data.")

        # Parse logs
        self.macs_data = dict()
        for f in self.find_log_files('macs2', filehandles=True):
            self.parse_macs(f)

        # Filter to strip out ignored sample names
        self.macs_data = self.ignore_samples(self.macs_data)

        if len(self.macs_data) == 0:
            log.debug("Could not find any MACS2 data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} logs".format(len(self.macs_data)))
        self.write_data_file(self.macs_data, 'multiqc_macs')

        # Add drop rate to the general stats table
        headers = OrderedDict()
        headers['treatment_redundant_rate'] = {
            'title': 'Treatment Redundancy',
            'description': 'Redundant rate in treatment',
            'max': 1,
            'min': 0,
            'format': '{:,.2f}',
            'scale': 'RdYlBu-rev'
        }
        headers['control_redundant_rate'] = {
            'title': 'Control Redundancy',
            'description': 'Redundant rate in control',
            'max': 1,
            'min': 0,
            'format': '{:,.2f}',
            'scale': 'RdYlBu-rev'
        }
        self.general_stats_addcols(self.macs_data, headers)

    def parse_macs(self, f):
        regexes = {
            'name': r"# name = (.+)$",
            'fragment_size': r"# fragment size is determined as (\d+) bps",
            'treatment_fragments_total': r"# total fragments in treatment: (\d+)",
            'treatment_fragments_after_filtering': r"# fragments after filtering in treatment: (\d+)",
            'treatment_max_duplicates': r"# maximum duplicate fragments in treatment = (\d+)",
            'treatment_redundant_rate': r"# Redundant rate in treatment: ([\d\.]+)",
            'control_fragments_total': r"# total fragments in control: (\d+)",
            'control_fragments_after_filtering': r"# fragments after filtering in control: (\d+)",
            'control_max_duplicates': r"# maximum duplicate fragments in control = (\d+)",
            'control_redundant_rate': r"# Redundant rate in control: ([\d\.]+)",
            'd': r"# d = (\d+)",
        }
        s_name = f['s_name']
        parsed_data = dict()
        for l in f['f']:
            for k, r in regexes.items():
                match = re.search(r, l)
                if match:
                    if k == 'name':
                        s_name = self.clean_s_name(match.group(1).strip(), f['root'])
                    else:
                        parsed_data[k] = float(match.group(1).strip())
            if not l.startswith('#') and l.strip():
                break
        if len(parsed_data) > 0:
            if s_name in self.macs_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.macs_data[s_name] = parsed_data

