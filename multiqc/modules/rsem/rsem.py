#!/usr/bin/env python

""" MultiQC module to parse output from Cutadapt """

from __future__ import print_function
import logging
import re
from distutils.version import StrictVersion
from collections import OrderedDict

from multiqc import config
#from multiqc.plots import linegraph
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    Cutadapt module class, parses stderr logs.
    Also understands logs saved by Trim Galore!
    (which contain cutadapt logs)
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Rsem', anchor='rsem',
        href='https://deweylab.github.io/RSEM/',
        info="RSEM (RNA-Seq by Expectation-Maximization) is a software package for"\
             "estimating gene and isoform expression levels from RNA-Seq data.")
    
         # Find and load any featureCounts reports
        self.rsem_mappable_data = dict()
        self.rsem_mapped_data = dict()
        self.rsem_keys  = list()
        
        for f in self.find_log_files('rsem'):
            
            if f['s_name'] in self.rsem_mapped_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            else :
                self.rsem_keys.append(f['s_name'])
            self.rsem_mapped_data[f['s_name']] = self.parse_rsem_report(f)
            
            self.add_data_source(f)
        # Filter to strip out ignored sample names
        self.rsem_mappable_data = self.ignore_samples(self.rsem_mappable_data)
        self.rsem_mapped_data = self.ignore_samples(self.rsem_mapped_data)

        if len(self.rsem_mapped_data) == 0  :
            log.info("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.rsem_mapped_data)))

        # Write parsed report data to a file
        self.write_data_file(self.rsem_mapped_data, 'multiqc_mappable_data')

        # Basic Stats Table
        self.rsem_stats_table()

        # Assignment bar plot
        self.add_section( plot = self.rsem_chart() )


    def parse_rsem_report (self, f):
        """ Parse the rsem cnt stat file. """
        data = dict()
        for l in f['f'].splitlines():
            s = l.split(" ")
            if len(s) > 3:
                #Line N0 N1 N2 N_tot   
                # N0, number of unalignable reads;
                # N1, number of alignable reads; 
                # N2, number of filtered reads due to too many alignments; N_tot = N0 + N1 + N2
                data['Unalignable'] = s[0]
                data['Alignable'] = s[1]
                data['Filtered'] = s[2]
                pass
            elif len(s) == 3:
                data['Unique'] = s[0]
                data['Multi'] = s[1]
                data['Uncertain'] = s[2]
            else:
                break
        return data
    
    def rsem_stats_table(self):
        """ Take the parsed stats from the rsem report and add them to the
        basic stats table at the top of the report """

        headers = OrderedDict()

        headers['Unalignable'] = {
            'title': '{} reads unalignable '.format(config.read_count_prefix),
            'description': 'number of unalignable reads ({})'.format(config.read_count_desc),
            'scale': 'PuBu',
            'modify': lambda x: float(x) * config.read_count_multiplier,
            'shared_key': 'read_count'
        }

        headers['Alignable'] = {
            'title': '{} reads alignable'.format(config.read_count_prefix),
            'description': 'number of alignable reads ({})'.format(config.read_count_desc),
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: float(x) * config.read_count_multiplier,
            'shared_key': 'read_count'
        }

        headers['Filtered'] = {
            'title': '{} reads filtered'.format(config.read_count_prefix),
            'description': 'number of filtered reads due to too many alignments ({})'.format(config.read_count_desc),
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: float(x) * config.read_count_multiplier,
            'shared_key': 'read_count'
        }

        headers['Unique'] = {
            'title': '{} reads aligned uniquely '.format(config.read_count_prefix),
            'description': 'number of reads aligned uniquely to a gene ({})'.format(config.read_count_desc),
            'scale': 'PuBu',
            'modify': lambda x: float(x) * config.read_count_multiplier,
            'shared_key': 'read_count'
        }
        headers['Multi'] = {
            'title': '{} reads aligned to multiple genes'.format(config.read_count_prefix),
            'description': 'number of reads aligned to multiple genes ({})'.format(config.read_count_desc),
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: float(x) * config.read_count_multiplier,
            'shared_key': 'read_count'
        }
        headers['Uncertain'] = {
            'title': '{} reads aligned to multiple locations'.format(config.read_count_prefix),
            'description': 'number of reads aligned  ({}) to multiple locations in the given reference sequences, which include isoform-level multi-mapping reads'.format(config.read_count_desc),
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: float(x) * config.read_count_multiplier,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.rsem_mapped_data, headers)


    def rsem_chart (self):
        """ Make the rsem assignment rates plot """
        
        # Config for the plot
        config = {
            'id': 'rsem_assignment_plot',
            'title': 'rsem Assignments',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        log.debug(self.rsem_mapped_data)
        log.debug(" " .join(self.rsem_keys))
        return bargraph.plot(self.rsem_mapped_data, ['Unique','Multi','Uncertain'],config)
        