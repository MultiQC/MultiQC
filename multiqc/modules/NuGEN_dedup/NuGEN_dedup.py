#!/usr/bin/env python

""" MultiQC module to parse output from FastQC
"""

from __future__ import print_function
#from __future__ import absolute_import
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
        super(MultiqcModule, self).__init__(name='NuGEN_dedup',
        anchor='NuGEN_dedup',
        href='sysbio.ucsd.edu',
        info="is a tool for deduplication of the reads")

        #find and load any featureCounts report
        self.dedup_data = dict()
        self.dedup_keys = list()

        for f in self.find_log_files('NuGEN_dedup'):
            fileName = f['s_name'][:-8]
            if fileName in self.dedup_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name'][:-8]))
            self.dedup_data[fileName] = self.parse_dedup_report(f['f'])
            self.dedup_keys = list(self.dedup_data[fileName].keys())
        
        if len(self.dedup_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.dedup_data)))
  
        #output stats file
        self.write_data_file(self.dedup_data, 'multiqc_dedup_cus')
        
        #add col to general stats table
        self.dedup_stats_table()
        
        #add section of bar_plot
        self.add_section (plot = self.dedup_frac_moltag_chart())
        self.add_section (plot = self.dedup_frac_position_chart()) 
    
    def parse_dedup_report (self,f):
        """ Parse the featureCounts log file. """
        dedup_keys = list()
        parsed_data = dict()
        infoList = []

        for l in f.splitlines():
            s = l.split()
            infoList.append(s)
        
        if len(infoList[0]) != len(infoList[1]):
            return "Incomplete table"
        
        for i in range (len(infoList[0])):
            if infoList[0][i] not in dedup_keys:
                dedup_keys.append(infoList[0][i])
            parsed_data[infoList[0][i]] = float(infoList[1][i])

        parsed_data['frac_moltag_dup'] = float(parsed_data['frac_moltag_dup'])*100
        parsed_data['moltag_unique_count'] = parsed_data['aligned_count']-parsed_data['moltag_dup_count']
        parsed_data['position_unique_count'] = parsed_data['aligned_count']-parsed_data['position_dup_count']
        return parsed_data

    def dedup_stats_table(self):
        headers = OrderedDict()
        headers['frac_moltag_dup'] = {
            'title': '% nudup_rate',
            'description':'NuGEN deduplication based on unique mapped reads',
            'max':100,
            'min':0,
            'suffix': '%',
            'scale': 'RdBu'
        }

        self.general_stats_addcols(self.dedup_data,headers)
    
    def dedup_frac_moltag_chart(self):
        keys_dedup = ['moltag_unique_count','moltag_dup_count']
        config = {
	    'id': 'dedup_frac_moltag_assignment_plot',
            'title': 'NuGEN_dedup_frac_moltag: Assignments',
            'ylab':'# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        return bargraph.plot(self.dedup_data, keys_dedup, config)

    def dedup_frac_position_chart(self):
        keys_dedup_position = ['position_unique_count','position_dup_count'] 
        config = {
	    'id': 'dedup_frac_position_assignment_plot',
	    'title': 'NuGEN_dedup_frac_position: Assignments',
	    'ylab':'# Reads',
            'cpswitch_counts_label': 'Number of Reads'
	}
        return bargraph.plot(self.dedup_data, keys_dedup_position, config)
