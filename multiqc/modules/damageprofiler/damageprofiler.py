#!/usr/bin/env python

""" MultiQC module to parse output from DamageProfiler """

from __future__ import print_function
from collections import OrderedDict
import logging
import json

from multiqc import config
from multiqc.plots import bargraph, linegraph
from multiqc.modules.base_module import BaseMultiqcModule

class MultiqcModule(BaseMultiqcModule):
    """
    damageprofiler module class
    """

    # Initialise the parent object
        super(MultiqcModule, self).__init__(name='damageprofiler',
        anchor='damageprofiler',
        href='https://github.com/Integrative-Transcriptomics/DamageProfiler',
        info="A Java based tool to determine damage patterns on ancient DNA as a replacement for mapDamage.")

        # Find and load 3pGtoAFreq Files
        self.3pGtoAfreq_data = dict() 
        for f in self.find_log_files('damageprofiler/threeprime'):
            self.parse3pG(f)
        
        # Find and load 5pCtoTFreq Files
        self.5pCtoTfreq_data = dict() 
        for f in self.find_log_files('damageprofiler/fiveprime'):
            self.parse5pC(f)
        
        # Find and load lgdist forward Files
        self.lgdist_fw_data = dict() 
        for f in self.find_log_files('damageprofiler/lgdistfw'):
            self.lgdistfw(f)

        # Find and load lgdist reverse Files
        self.lgdist_rv_data = dict() 
        for f in self.find_log_files('damageprofiler/lgdistrv'):
            self.lgdistrv(f)

        # Filter to strip out ignored sample names
        self.3pGtoAfreq_data         =   self.ignore_samples(self.3pGtoAfreq_data)
        self.5pCtoTfreq_data          =   self.ignore_samples(self.5pCtoTfreq_data)
        self.lgdist_fw_data   =   self.ignore_samples(self.lgdist_fw_data)
        self.lgdist_rv_data      =   self.ignore_samples(self.lgdist_rv_data)

        # Warning when no files are found
        if max(len(self.3pGtoAfreq_data), len(self.5pCtoTfreq_data), len(self.lgdist_fw_data), len(self.lgdist_rv_data)) == 0:
            raise UserWarning
        
        
        
        