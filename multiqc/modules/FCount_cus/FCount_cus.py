#/usr/bin/env python

""" MultiQC module to parse output from FastQC
"""
from __future__ import print_function
from collections import OrderedDict
import logging
import re
import csv

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(name='FCount_cus',
        anchor='FCount_cus',
        href='sysbio.ucsd.edu',
        info="is a tool for count the number of genes that appear specific times")

        self.FCount_data = dict()
        self.FCount_keys = list()
        self.FCount_report = list()

        #search for specific files 
        for f in self.find_log_files('FCount_cus'):
            fileName = f['s_name']
            
            if fileName in self.FCount_report:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))

            self.FCount_report.append(fileName)
            self.FCount_data = self.parse_FCount_report(f['f'])
            self.FCount_keys = self.FCount_data.keys()
   
        #raise user warinings and information   
        if len(self.FCount_keys) == 0:
            raise UserWarning
        if len(self.FCount_report) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.FCount_report)))

        #output stats file
        self.write_data_file(self.FCount_data, 'multiqc_geneCount_cus')

        #add col to general stats table
        self.FCount_stats_table()

        #add report section
        self.add_section(plot = self.FCount_barPlot_0())

    def parse_FCount_report(self,f):
        #Parse the mergeCount table
        parsed_data = dict()
        infoList = []
        count = 0
        for l in f.splitlines():
            s = l.split()
            fileName = s[0]
            
            #skip the firt line (title)
            if fileName =='sample':
                continue
            
            count_info = dict()
            count_info["counts > 0"] = int(s[1])
            count_info["counts > 10"] = int(s[2])
            count_info["counts > 100"] = int(s[3])
            
            #information for report section
            count_info["counts:0-10"] = count_info["counts > 0"]-count_info["counts > 10"]-count_info["counts > 100"]
            count_info["counts:10-100"] = count_info["counts > 10"]-count_info["counts > 100"]
            
            parsed_data[fileName] = count_info
        
        return parsed_data

    def FCount_stats_table(self):
        headers = OrderedDict()
        headers['counts > 0'] = {
            'title':'count > 0',
            'description':'Count the genes that appear one or more times after the NuGEN deduplication process',
            'scale': 'RdYlGn-rev',
            'placement' : 1.0,
        }
        headers['counts > 10'] = {
            'title':'count > 10',
            'description':'Count the genes that appear ten or more times after the NuGEN deduplication process',
            'scale': 'RdBu',
            'hidden': True,
            'placement' : 2.0,
        }
        headers['counts > 100'] = {
            'title':'count > 100',
            'description':'Count the genes that appear one hundred or more times after the NuGEN deduplication process',
            'scale': 'RdBu',
            'hidden': True,
            'placement' : 3.0,
        }
        self.general_stats_addcols(self.FCount_data,headers)

    def FCount_barPlot_0(self):
        keys_FCount = ['counts:0-10','counts:10-100','counts > 100']
        config = {
            'id': 'Count_Distribution_StatsPlot',
            'title': 'Counts_Distribution: StatsPlot',
            'ylab':'# Counts',
            'cpswitch_counts_label': 'Number of Genes'
    }
        return bargraph.plot(self.FCount_data, keys_FCount, config)
