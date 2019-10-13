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
        self.sexdet_nraut = dict()
        self.sexdet_nrX = dict()
        self.sexdet_nrY = dict()
        self.sexdet_rateErrX = dict()
        self.sexdet_rateErrY = dict()
        self.sexdet_rateX = dict()
        self.sexdet_rateY = dict() 
        self.sexdet_snpsauto = dict()
        self.sexdet_Xsnps = dict()
        self.sexdet_Ysnps = dict()

        # Find and load JSON file
        for f in self.find_log_files('sexdeterrmine',filehandles=True):
            self.parseJSON(f)
        
        #Filter to strip out ignored sample names
        #self.sexdet_nraut = self.ignore_samples(self.sexdet_nraut)
        #self.sexdet_nrX = self.ignore_samples(self.sexdet_nrX)
        #self.sexdet_nrY = self.ignore_samples(self.sexdet_nrY)
        #self.sexdet_rateErrX = self.ignore_samples(self.sexdet_rateErrX)
        #self.sexdet_rateErrY = self.ignore_samples(self.sexdet_rateErrY)
        #self.sexdet_rateX = self.ignore_samples(self.sexdet_rateX)
        #self.sexdet_rateY = self.ignore_samples(self.sexdet_rateY)
        #self.sexdet_snpsauto = self.ignore_samples(self.sexdet_snpsauto)
        #self.sexdet_Xsnps = self.ignore_samples(self.sexdet_Xsnps)
        #self.sexdet_Ysnps = self.ignore_samples(self.sexdet_Ysnps)

        self.write_data_file(self.sexdet_nraut, 'multiqc_sexdeter_nraut_metrics')


        #Parse our nice little JSON file
    def parseJSON(self, f):
        """ Parse the JSON output from SexDeterrmine and save the summary statistics """
        try:
            parsed_json = json.load(f['f'])
            pp = pprint.PrettyPrinter(indent=4)
            pp.pprint(parsed_json)
        except Exception as e:
            print(e)
            log.warn("Could not parse SexDeterrmine JSON: '{}'".format(f['fn']))
            return None

        pp.pprint(parsed_json['/Users/alexanderpeltzer/Downloads/BOO002.sorted.bam'])
        #Get sample name from JSON first
        for k in parsed_json:
            if (k != 'Metadata'):
                try:  
                    s_name = self.clean_s_name(k,'')
                    self.add_data_source(f, s_name)
                    self.sexdet_nraut[s_name] = parsed_json[k]['NR Aut']
                    self.sexdet_nrX[s_name] = parsed_json[k]['NrX']
                    self.sexdet_nrY[s_name] = parsed_json[k]['NrY']
                    self.sexdet_rateErrX[s_name] = float(parsed_json[k]['RateErrX'])
                    self.sexdet_rateErrY[s_name] = float(parsed_json[k]['RateErrY'])
                    self.sexdet_rateX[s_name] = float(parsed_json[k]['RateX'])
                    self.sexdet_rateY[s_name] = float(parsed_json[k]['RateY'])
                    self.sexdet_snpsauto[s_name] = parsed_json[k]['Snps Autosomal']
                    self.sexdet_Xsnps[s_name] = parsed_json[k]['XSnps']
                    self.sexdet_Ysnps[s_name] = parsed_json[k]['YSnps']
                except ValueError as error:
                    print("Something isn't right with this error:", error)
            else:
                continue