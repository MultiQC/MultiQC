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
        super(MultiqcModule, self).__init__(name='DamageProfiler:',
        anchor='damageprofiler',
        href='https://github.com/Integrative-Transcriptomics/DamageProfiler',
        info="a tool to determine damage patterns on ancient DNA.")

        # Find and load 3pGtoAFreq Files
        self.threepGtoAfreq_data = dict() 
        for f in self.find_log_files('damageprofiler/threeprime'):
            self.parse3pG(f)

        # Find and load 5pCtoTFreq Files
        self.fivepCtoTfreq_data = dict() 
        for f in self.find_log_files('damageprofiler/fiveprime'):
            self.parse5pC(f)

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
        self.fivepCtoTfreq_data          =   self.ignore_samples(self.fivepCtoTfreq_data)
        #self.lgdist_fw_data   =   self.ignore_samples(self.lgdist_fw_data)
        #self.lgdist_rv_data      =   self.ignore_samples(self.lgdist_rv_data)

        # Warning when no files are found
        #if max(len(self.threepGtoAfreq_data), len(self.5pCtoTfreq_data), len(self.lgdist_fw_data), len(self.lgdist_rv_data)) == 0:
        #    raise UserWarning

        # Write parsed data to a file
        self.write_data_file(self.threepGtoAfreq_data, 'multiqc_damageprofiler_3pGtoAfreq')
        self.write_data_file(self.fivepCtoTfreq_data, 'multiqc_damageprofiler_5pCtoTfreq')
        #self.write_data_file(self.lgdist_fw_data, 'multiqc_damageprofiler_lgdist_fw')
        #self.write_data_file(self.lgdist_rv_data, 'multiqc_damageprofiler_lgdist_rv')

        # Basic Stats Table
        self.damageprofiler_general_stats_table_3p()
        #self.general_stats_addcols(self.threepGtoAfreq_data)

        # Add plots
        if len(self.threepGtoAfreq_data) > 0:
            self.add_section ( 
                name = '3\' Misincorporation Plot',
                plot = self.threeprime_plot()
            )
        if len(self.fivepCtoTfreq_data) > 0:
            self.add_section (
                name = '5\' Misincorporation Plot',
                plot = self.fiveprime_plot()
            )


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
            # Create tuples out of entries
            pos = list(range(1,len(parsed_data.get(key))))
            tuples = list(zip(pos,parsed_data.get(key)))
            # Get a dictionary out of it
            data = dict((x, y) for x, y in tuples)
            self.threepGtoAfreq_data[s_name] = data
        else: 
            log.debug('No valid data {} in 3pGtoA report'.format(f['fn']))
            return None
    
    #Parse a 5pCtoTfreq file
    def parse5pC(self, f):

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

            if s_name in self.fivepCtoTfreq_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            # Create tuples out of entries
            pos = list(range(1,len(parsed_data.get(key))))
            tuples = list(zip(pos,parsed_data.get(key)))
            # Get a dictionary out of it
            data = dict((x, y) for x, y in tuples)
            self.fivepCtoTfreq_data[s_name] = data
        else: 
            log.debug('No valid data {} in 5pCtoT report'.format(f['fn']))
            return None
    
    #### Tables from here on 
    def damageprofiler_general_stats_table_3p(self):
        """ Take the parsed stats from the DamageProfiler and add it to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['damageprofiler'] = {
            'title': '3p G to A',
            'description': 'Percentage of misincorporated G to A substitutions.',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'OrRd',
            'format': '{:,.0f}',
            'modify': lambda x: x * 100.0
        }
        self.general_stats_addcols(self.threepGtoAfreq_data, headers)


    #### Plotting from here on


    #Linegraph plot for 3pGtoA
    def threeprime_plot(self):
        """Generate a 3' GtoA linegraph plot"""

        data = dict()
        for s_name in self.threepGtoAfreq_data:
            try:
                data[s_name] = {int(d): float (self.threepGtoAfreq_data[s_name][d])*100 for d in self.threepGtoAfreq_data[s_name]}
            except KeyError:
                pass
        if len(data) == 0:
            log.debug('No valid data for 3\' G to A input!')
            return None

        config = {
            'id': 'threeprime_misinc_plot',
            'title': 'DamageProfiler: 3P G to A Misincorporation plot',
            'ylab': '% G to A substituted',
            'xlab': 'Nucleotide Position from 3\'',
            'ymin': 0,
            'xmin': 1
        }

        return linegraph.plot(data,config)

    #Linegraph plot for 3pGtoA
    def fiveprime_plot(self):
        """Generate a 5' CtoT linegraph plot"""

        data = dict()
        for s_name in self.fivepCtoTfreq_data:
            try:
                data[s_name] = {int(d): float (self.fivepCtoTfreq_data[s_name][d])*100 for d in self.fivepCtoTfreq_data[s_name]}
            except KeyError:
                pass
        if len(data) == 0:
            log.debug('No valid data for 5\' C to T input!')
            return None

        config = {
            'id': 'fiveprime_misinc_plot',
            'title': 'DamageProfiler: 5P C to T Misincorporation plot',
            'ylab': '% C to T substituted',
            'xlab': 'Nucleotide Position from 5\'',
            'ymin': 0,
            'xmin': 1
        }

        return linegraph.plot(data,config)

        
        