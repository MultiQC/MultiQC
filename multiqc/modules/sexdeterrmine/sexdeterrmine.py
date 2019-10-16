#!/usr/bin/env python

""" MultiQC module to parse output from SexdetErrmine """

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import re
import json
import pprint

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ SexDeterrmine module """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='SexDetErrmine', anchor='sexdeterrmine',
        href="https://github.com/TCLamnidis/Sex.DetERRmine",
        info="A python script to calculate the relative coverage of X and Y chromosomes, and their associated error bars, from the depth of coverage at specified SNPs. ")

        # Find and load any DeDup reports
        self.sexdet_data = dict()
        

        # Find and load JSON file
        for f in self.find_log_files('sexdeterrmine',filehandles=True):
            self.parseJSON(f)
        
        # Empty dictionary for sending to table
        dict_to_plot = OrderedDict()

        pp = pprint.PrettyPrinter(indent=4) 
        pp.pprint(self.sexdet_data)
        for k in self.sexdet_data:
            if (k != 'Metadata'):
                try:  
                    s_name = self.clean_s_name(k,f['root'])
                    self.add_data_source(f, s_name)
                    dict_to_plot[s_name]['NR Aut'] = self.sexdet_data[k]['NR Aut']
                    dict_to_plot[s_name]['NrX'] = self.sexdet_data[k]['NrX']
                    dict_to_plot[s_name]['NrY'] = self.sexdet_data[k]['NrY']
                    dict_to_plot[s_name]['RateErrX'] = self.sexdet_data[k]['RateErrX']
                    dict_to_plot[s_name]['RateErrY'] = self.sexdet_data[k]['RateErrY']
                    dict_to_plot[s_name]['RateX'] = self.sexdet_data[k]['RateX']
                    dict_to_plot[s_name]['RateY'] = self.sexdet_data[k]['RateY']
                    dict_to_plot[s_name]['Snps Autosomal'] = self.sexdet_data[k]['Snps Autosomal']
                    dict_to_plot[s_name]['XSnps'] = self.sexdet_data[k]['XSnps']
                    dict_to_plot[s_name]['YSnps'] = self.sexdet_data[k]['YSnps']
                    self.addSummaryMetrics(dict_to_plot)
                    dict_to_plot = OrderedDict()
                except ValueError as error:
                    print("Something isn't right with this error:", error)
            else:
                continue

       
        
        
        self.write_data_file(self.sexdet_data, 'multiqc_sexdeter_metrics')


        #Parse our nice little JSON file
    def parseJSON(self, f):

        """ Parse the JSON output from SexDeterrmine and save the summary statistics """
        try:
            parsed_json = json.load(f['f'])
            self.sexdet_data = parsed_json
        except Exception as e:
            print(e)
            log.warn("Could not parse SexDeterrmine JSON: '{}'".format(f['fn']))
            return None

    def addSummaryMetrics(self, dict_to_plot):
        """ Take the parsed stats from SexDetErrmine and add it to the main plot """

        headers = OrderedDict()
        headers['nraut'] = {
            'title': '# Autosomal Pos',
            'description': 'The number of reads covering positions on the autosome.',
            'scale': 'PuBu',
            'format': '{:,.2f}',
            'hidden': True
        }
        headers['nrX'] = {
            'title': '# X Pos',
            'description': 'The number of reads covering positions on Chromosome X.',
            'scale': 'YlGnBu',
            'format': '{:,.2f}',
            'hidden': True
        }
        headers['nrY'] = {
            'title': '# Y Pos',
            'description': 'The number of reads covering positions on Chromosome Y.',
            'scale': 'YlGnBu',
            'format': '{:,.2f}',
            'hidden': True
        }
        headers['rateErrX'] = {
            'title': 'rateErrX',
            'description': 'Rate of Error for Chr X',
            'scale': 'PuBuGn',
            'format': '{:,.2f}',
            'hidden': True
        }
        headers['rateErrY'] = {
            'title': 'rateErrY',
            'description': 'Rate of Error for Chr Y',
            'scale': 'PuBuGn',
            'format': '{:,.2f}',
            'hidden': True
        }
        headers['rateX'] = {
            'title': 'rateX',
            'description': 'Number of positions on Chromosome X vs Autosomal positions.',
            'scale': 'PuBuGn',
            'format': '{:,.2f}',
            'hidden': True
        }
        headers['rateY'] = {
            'title': 'rateY',
            'description': 'Number of positions on Chromosome Y vs Autosomal positions.',
            'scale': 'PuBuGn',
            'format': '{:,.2f}',
            'hidden': True
        }
        headers['snpsauto'] = {
            'title': 'Pos on Auto',
            'description': 'Total number of autosomal positions. When supplied with a BED file, this includes only positions specified there.',
            'scale': 'PuBuGn',
            'format': '{:,.2f}',
            'hidden': True
        } 
        headers['XSnps'] = {
            'title': 'Pos on X',
            'description': 'Total number of positions on Chromosome X. When supplied with a BED file, this includes only positions specified there.',
            'scale': 'PuBuGn',
            'format': '{:,.2f}',
            'hidden': True
        } 
        headers['YSnps'] = {
            'title': 'Pos on Y',
            'description': 'Total number of positions on Chromosome Y. When supplied with a BED file, this includes only positions specified there.',
            'scale': 'PuBuGn',
            'format': '{:,.2f}',
            'hidden': True
        } 

        self.general_stats_addcols(dict_to_plot, headers)
