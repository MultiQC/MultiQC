#!/usr/bin/env python


""" MultiQC module to parse output from sickle """

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
import logging
import re
from multiqc.plots import bargraph
from collections import OrderedDict

log = logging.getLogger(__name__)
class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(name='sickle', anchor='sickle',
                                            href="https://github.com/najoshi/sickle",
                                            info="A windowed adaptive trimming tool for FASTQ files using quality.")

        self.sickle_data = dict()
        #parse list of log files
        for f in self.find_log_files('sickle'):
            self.sickle_data[f['s_name']] = self.parse_logs(f['f'])
            self.add_data_source(f)
        #no file found
        if len(self.sickle_data) == 0:
            raise UserWarning
        self.sickle_data = self.ignore_samples(self.sickle_data)
        self.sickle_general_stats_table()
        self.write_data_file(self.sickle_data, 'multiqc_sickle')
        self.read_count_plot()

    def parse_logs(self, f):
        """ Parse the Sickle standard output """
        regexes = {
                'reads_paired_kept': "FastQ paired records kept: ([\d,]+) .*",
                'reads_single_kept': "FastQ single records kept: ([\d,]+).*",
                'reads_paired_discarded': "FastQ paired records discarded: ([\d,]+) .*",
                'reads_single_discarded': "FastQ single records discarded: ([\d,]+) .*",
                'reads_kept': "FastQ records kept: ([\d,]+)",
                'reads_discarded': "FastQ records discarded: ([\d,]+)"
            }
        data = {}
        for l in f.splitlines():
            s = l.split(':')
            # Search regexes for overview stats# Search regexes for overview stats
            for k, r in regexes.items():
                match = re.search(r, l)
                if match:
                    data[k] = int(match.group(1).replace(',', ''))
        #if paired data : compute total for general stats table
        if "reads_kept" not in data.keys() :
            data["reads_kept"] = data["reads_paired_kept"] + data["reads_single_kept"]

        if "reads_discarded" not in data.keys() :
            data["reads_discarded"] = data["reads_paired_discarded"] + data["reads_single_discarded"]

        return data

    def sickle_general_stats_table(self):
        """ Take the parsed stats from the sickle report and add
        number of kept reads into the generam stats table """

        headers = {}
        headers['reads_kept'] = {
            'title': '{} kept reads (sickle)'.format(config.read_count_prefix),
            'description': 'Number of reads kept after sikle ({})'.format(config.read_count_desc),
            'modify': lambda x: x * config.read_count_multiplier
        }
        self.general_stats_addcols(self.sickle_data, headers)

    def read_count_plot (self):
        """ Stacked bar plot showing counts of reads """
        pconfig = {
            'id': 'sickle_sequence_counts_plot',
            'title': 'Sickle: Sequence Counts',
            'ylab': 'Number of reads',
            'cpswitch_counts_label': 'Number of reads',
            'hide_zero_cats': False
        }
        pdata = dict()
        has_paired = False
        for s_name in self.sickle_data:
            pd = self.sickle_data[s_name]
            pdata[s_name] = dict()

            if "reads_paired_kept" in pd.keys() : # paired library 4 values to add in barplot
                pdata[s_name]['reads_paired_kept'] = pd['reads_paired_kept']
                pdata[s_name]['reads_single_kept'] = pd['reads_single_kept']
                pdata[s_name]['reads_paired_discarded'] =  pd['reads_paired_discarded']
                pdata[s_name]['reads_single_discarded'] =  pd['reads_single_discarded']
                has_paired = True
            else : # single lib only 2 values
                pdata[s_name]['reads_single_kept'] = pd['reads_kept']
                pdata[s_name]['reads_single_discarded'] =  pd['reads_discarded']

        pcats = OrderedDict()
        if has_paired :
            pcats['reads_paired_kept']={
                    'name': 'Paired reads kept',
                    'color': '#1f78b4'
            }
        pcats['reads_single_kept']={
                'name': 'Single reads kept',
                'color': '#a6cee3'
            }

        if has_paired :
            pcats['reads_paired_discarded']= {
                'name': 'Paired reads discarded',
                'color': '#ff7f00'
            }
        pcats['reads_single_discarded']={
                'name': 'Single reads discarded',
                'color': '#fdae61'
            }

        self.add_section (
            plot = bargraph.plot(pdata, pcats, pconfig)
        )
