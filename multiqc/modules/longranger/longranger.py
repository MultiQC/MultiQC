#!/usr/bin/env python

""" MultiQC module to parse output from Longranger """

from __future__ import print_function
from collections import OrderedDict
import logging
import re
import os
import csv

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ Longranger module """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Longranger', anchor='longranger',
        href="https://www.10xgenomics.com/",
        info="A set of analysis pipelines that perform sample demultiplexing, "
        "barcode processing, alignment, quality control, variant calling, phasing, "
        "and structural variant calling.")

        self.longranger_data = dict()
        self.paths_dict = dict()
        for f in self.find_log_files('longranger/invocation'):
            sid = self.parse_invocation(f['f'])
            self.paths_dict[os.path.basename(f['root'])] = sid
        
        running_name = 1
        for f in self.find_log_files('longranger/summary'):
            data = self.parse_summary(f['f'])
            updir, _ = os.path.split(f['root'])
            base_updir = os.path.basename(updir)
            sid = "Longranger_{}".format(running_name)
            if base_updir in self.paths_dict.keys():
                sid = self.paths_dict[base_updir]
            else:
                running_name += 1
            
            self.longranger_data[sid] = data 

        # Filter to strip out ignored sample names
        self.longranger_data = self.ignore_samples(self.longranger_data)

        if len(self.longranger_data) == 0:
            raise UserWarning

        self.headers = OrderedDict()
        self.headers['molecule_length_mean'] = {
                'description': 'placeholder'
        }
        self.headers['longranger_version'] = {
                'description': 'placeholder'
        }
        self.headers['mean_depth'] = {
                'description': 'placeholder'
        }
        self.headers['mapped_reads'] = {
                'description': 'placeholder'
        }
        self.headers['corrected_loaded_mass_ng'] = {
                'description': 'placeholder'
        }
        self.headers['snps_phased'] = {
                'description': 'placeholder'
        }
        self.headers['n50_phase_block'] = {
                'description': 'placeholder'
        }



        ### Write the report
        self.write_data_file(self.longranger_data, 'multiqc_longranger')
        config_table = {
            'id': 'longranger_table',
            'namespace': 'longranger'
        }
        self.add_section (
            name = 'Longranger statistics',
            anchor = 'longranger-table',
            description = 'placeholder',
            helptext = 'placeholder',
            plot = table.plot(self.longranger_data, self.headers, config_table)
        )


        log.info("Found {} reports".format(len(self.longranger_data)))

        # Write parsed report data to a file
        self.write_data_file(self.longranger_data, 'multiqc_longranger')


    def parse_invocation(self, content):
        sid_pat = re.compile('    sample_id = \"(.*)\",')
        
        sid = None
        for l in content.splitlines():
            sid_m = re.match(sid_pat,l)
            if sid_m is not None:
                sid = sid_m.groups()[0]
        return sid

    def parse_summary(self, content):
        
        out_dict = OrderedDict()
        lines = content.splitlines()
        data = list(zip(lines[0].strip().split(','), lines[1].strip().split(',')))
        for i,j in data:
            out_dict[i] = j
            
        return out_dict
