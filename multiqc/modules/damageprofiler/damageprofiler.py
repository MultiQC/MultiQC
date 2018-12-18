#!/usr/bin/env python

""" MultiQC module to parse output from DamageProfiler """

from __future__ import print_function
from collections import OrderedDict
import logging
import yaml
import pprint

from multiqc import config
from multiqc.plots import  linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):
        """
        damageprofiler module class
        """

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='damageprofiler',
        anchor='damageprofiler',
        href='https://github.com/Integrative-Transcriptomics/DamageProfiler',
        info="A Java based tool to determine damage patterns on ancient DNA as a replacement for mapDamage.")

        # Find and load 3pGtoAFreq Files
        self.threepGtoAfreq_data = dict() 
        for f in self.find_log_files('damageprofiler/threeprime'):
            self.parse3pG(f)

        # Find and load 5pCtoTFreq Files
        #self.5pCtoTfreq_data = dict() 
        #for f in self.find_log_files('damageprofiler/fiveprime'):
        #    self.parse5pC(f)

        # Find and load lgdist forward Files
        #self.lgdist_fw_data = dict() 
        #for f in self.find_log_files('damageprofiler/lgdistfw'):
        #    self.lgdistfw(f)

        # Find and load lgdist reverse Files
        #self.lgdist_rv_data = dict() 
        #for f in self.find_log_files('damageprofiler/lgdistrv'):
        #    self.lgdistrv(f)

        # Filter to strip out ignored sample names
        self.threepGtoAfreq_data         =   self.ignore_samples(self.threepGtoAfreq_data)
        #self.5pCtoTfreq_data          =   self.ignore_samples(self.5pCtoTfreq_data)
        #self.lgdist_fw_data   =   self.ignore_samples(self.lgdist_fw_data)
        #self.lgdist_rv_data      =   self.ignore_samples(self.lgdist_rv_data)

        # Warning when no files are found
        #if max(len(self.threepGtoAfreq_data), len(self.5pCtoTfreq_data), len(self.lgdist_fw_data), len(self.lgdist_rv_data)) == 0:
        #    raise UserWarning

        # Write parsed data to a file
        #self.write_data_file(self.threepGtoAfreq_data, 'multiqc_damageprofiler_3pGtoAfreq')
        #self.write_data_file(self.5pCtoTfreq_data, 'multiqc_damageprofiler_5pCtoTfreq')
        #self.write_data_file(self.lgdist_fw_data, 'multiqc_damageprofiler_lgdist_fw')
        #self.write_data_file(self.lgdist_rv_data, 'multiqc_damageprofiler_lgdist_rv')
    
    #Parse a 3pGtoAfreq file
    def parse3pG(self, f):

        try:
            # Parsing as OrderedDict is slightly messier with YAML
            # http://stackoverflow.com/a/21048064/713980
            # Copied over from custom_content.py - thanks @ewels!
            def dict_constructor(loader, node):
                return OrderedDict(loader.construct_pairs(node))
            yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, dict_constructor)
            parsed_data = yaml.load(f['f'])
        except Exception as e:
            log.warning("Error parsing YAML file '{}' (probably invalid YAML)".format(f['fn']))
            log.warning("YAML error: {}".format(e))
            return None
        
        #Sample name is always the only key in the YAML
        if parsed_data is not None:
            key = list(parsed_data.keys())[0] 
            s_name = self.clean_s_name(list(parsed_data.keys())[0],'')

            if s_name in self.threepGtoAfreq_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
        
            self.threepGtoAfreq_data[s_name] = parsed_data.get(key)
            pprint.pprint(self.threepGtoAfreq_data)
        else: 
            log.debug('No valid data {} in 3pGtoA report'.format(f['fn']))
            return None


        
        