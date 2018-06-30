#!/usr/bin/env python

""" MultiQC module to parse output from MACS2 """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import bargraph
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
            raise UserWarning

        log.info("Found {} logs".format(len(self.macs_data)))
        self.write_data_file(self.macs_data, 'multiqc_macs')

        self.macs_general_stats()
        self.macs_filtered_reads_plot()


    def parse_macs(self, f):
        regexes = {
            'name': r"# name = (.+)$",
            'fragment_size': r"# (?:fragment|tag) size is determined as (\d+) bps",
            'treatment_fragments_total': r"# total (?:fragments|tags) in treatment: (\d+)",
            'treatment_fragments_after_filtering': r"# (?:fragments|tags) after filtering in treatment: (\d+)",
            'treatment_max_duplicates': r"# maximum duplicate (?:fragments|tags at the same position) in treatment = (\d+)",
            'treatment_redundant_rate': r"# Redundant rate in treatment: ([\d\.]+)",
            'control_fragments_total': r"# total (?:fragments|tags) in control: (\d+)",
            'control_fragments_after_filtering': r"# (?:fragments|tags) after filtering in control: (\d+)",
            'control_max_duplicates': r"# maximum duplicate (?:fragments|tags at the same position) in control = (\d+)",
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

    def macs_general_stats(self):
        """ Add columns to General Statistics table """
        headers = OrderedDict()
        headers['d'] = {
            'title': 'Fragment Length',
            'min': 0,
            'format': '{:,.0f}'
        }
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

    def macs_filtered_reads_plot(self):
        """ Plot of filtered reads for control and treatment samples """
        data = dict()
        req_cats = ['control_fragments_total', 'control_fragments_after_filtering', 'treatment_fragments_total', 'treatment_fragments_after_filtering']
        for s_name, d in self.macs_data.items():
            if all([c in d for c in req_cats]):
                data['{}: Control'.format(s_name)] = dict()
                data['{}: Treatment'.format(s_name)] = dict()
                data['{}: Control'.format(s_name)]['fragments_filtered'] = d['control_fragments_total'] - d['control_fragments_after_filtering']
                data['{}: Control'.format(s_name)]['fragments_not_filtered'] = d['control_fragments_after_filtering']
                data['{}: Treatment'.format(s_name)]['fragments_filtered'] = d['treatment_fragments_total'] - d['treatment_fragments_after_filtering']
                data['{}: Treatment'.format(s_name)]['fragments_not_filtered'] = d['treatment_fragments_after_filtering']

        # Check that we have something to plot
        if len(data) == 0:
            return

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['fragments_not_filtered'] = { 'color': '#437BB1', 'name': 'Remaining fragments' }
        keys['fragments_filtered'] =     { 'color': '#B1084C', 'name': 'Filtered fragments' }

        # Config for the plot
        pconfig = {
            'id': 'macs2_filtered',
            'title': 'MACS2: Filtered Fragments',
            'ylab': '# Fragments',
            'cpswitch_counts_label': 'Number of Fragments',
            'hide_zero_cats': False
        }

        self.add_section(
            plot = bargraph.plot(data, keys, pconfig)
        )
