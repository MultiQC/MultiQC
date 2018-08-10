#!/usr/bin/env python

""" MultiQC module to parse Stacks denovo output"""

from __future__ import print_function
from collections import OrderedDict
import logging
import re
import os
from multiqc import config
from multiqc.plots import table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):
        super(MultiqcModule, self).__init__(name='Stacks', anchor='stacks',
        href="http://catchenlab.life.illinois.edu/stacks/",
        info="A software for analyzing restriction enzyme-based data (e.g. RAD-seq).")

        self.gstacks_data = OrderedDict()
        self.gsheaders = OrderedDict()
        self.gsheaders['n_loci'] = {
                'title': '# loci'
        }
        self.gsheaders['n_used_fw_reads'] = {
                'title': '# reads'
        }
        self.gsheaders['mean_cov'] = {
                'title': 'cov'
        }
        self.gsheaders['mean_cov_ns'] = {
                'title': 'weighted cov'
        }
        # Parse gstacks data
        for f in self.find_log_files('stacks/gstacks'):
            s_name = self.clean_s_name(os.path.basename(f['root']), os.path.dirname(f['root']))
            # todo: append data / change names for multiple runs
            self.gstacks_data = self.parse_gstacks(f['f'], s_name, f)

        # Parse populations data
        for f in self.find_log_files('stacks/populations'):
            s_name = self.clean_s_name(os.path.basename(f['root']), os.path.dirname(f['root']))



        # Write parsed report data to a file
        self.write_data_file(self.gstacks_data, 'multiqc_stacks')

        ### Write the table
        config_table = {
            'id': 'gstacks_table',
            'namespace': 'stacks'
        }
        self.add_section (
            name = 'gstacks',
            anchor = 'stacks-gstacks',
            description = 'placeholder',
            helptext = '''placeholder''',
            plot = table.plot(self.gstacks_data, self.gsheaders, config_table)
        )


    def parse_gstacks(self, file_contents, s_name=None, f=None):

        headers = None
        content = None
        out_dict = OrderedDict()

        for l in file_contents.splitlines():
            if l.startswith("sample	n_loci	n_used_fw_reads"):
                headers = ["n_loci", "n_used_fw_reads", "mean_cov","mean_cov_ns"]
            elif l.startswith("END effective_coverages_per_sample"):
                break
            elif headers is not None:
                cdict = OrderedDict().fromkeys(headers)
                content = l.split("\t")
                for i in range(1, len(content)):
                    cdict[list(cdict.keys())[i-1]] = content[i]
                out_dict[content[0]] = cdict


        return out_dict

    def parse_populations(self, file_contents, s_name=None, f=None):

        out_dict = OrderedDict()
        
