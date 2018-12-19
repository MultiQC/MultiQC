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
            self.parseFreqPlot(f,self.threepGtoAfreq_data)

        # Find and load 5pCtoTFreq Files
        self.fivepCtoTfreq_data = dict() 
        for f in self.find_log_files('damageprofiler/fiveprime'):
            self.parseFreqPlot(f,self.fivepCtoTfreq_data)

        # Find and load lgdist forward Files
        self.lgdist_fw_data = dict() 
        for f in self.find_log_files('damageprofiler/lgdistfw'):
            self.parselgDist(f,self.lgdist_fw_data)

        # Find and load lgdist reverse Files
        self.lgdist_rv_data = dict() 
        for f in self.find_log_files('damageprofiler/lgdistrv'):
            self.parselgDist(f,self.lgdist_rv_data)

        # Filter to strip out ignored sample names
        self.threepGtoAfreq_data         =   self.ignore_samples(self.threepGtoAfreq_data)
        self.fivepCtoTfreq_data          =   self.ignore_samples(self.fivepCtoTfreq_data)
        self.lgdist_fw_data   =   self.ignore_samples(self.lgdist_fw_data)
        self.lgdist_rv_data      =   self.ignore_samples(self.lgdist_rv_data)

        # Warning when no files are found
        #if max(len(self.threepGtoAfreq_data), len(self.5pCtoTfreq_data), len(self.lgdist_fw_data), len(self.lgdist_rv_data)) == 0:
        #    raise UserWarning

        # Write parsed data to a file
        self.write_data_file(self.threepGtoAfreq_data, 'multiqc_damageprofiler_3pGtoAfreq')
        self.write_data_file(self.fivepCtoTfreq_data, 'multiqc_damageprofiler_5pCtoTfreq')
        self.write_data_file(self.lgdist_fw_data, 'multiqc_damageprofiler_lgdist_fw')
        self.write_data_file(self.lgdist_rv_data, 'multiqc_damageprofiler_lgdist_rv')

        # Basic Stats Table, use generic function to add data to general table
        self.dmgprof_misinc_stats(self.threepGtoAfreq_data, '3p G to A', 'Percentage of misincorporated G to A substitutions.')
        self.dmgprof_misinc_stats(self.fivepCtoTfreq_data, '5P C to T', 'Percentage of misincorporated C to T substitutions.')

        # Add plots
        if len(self.threepGtoAfreq_data) > 0:
            self.add_section ( 
                name = '3P Misincorporation Plot',
                plot = self.threeprime_plot()
            )
        if len(self.fivepCtoTfreq_data) > 0:
            self.add_section (
                name = '5P Misincorporation Plot',
                plot = self.fiveprime_plot()
            )
        if len(self.lgdist_fw_data) > 0 and len(self.lgdist_rv_data) > 0:
            self.add_section (
                name = 'Forward read length distribution',
                plot = self.lgdistplot(self.lgdist_fw_data, 'Forward')
            )
            self.add_section (
                name = 'Reverse read length distribution',
                plot = self.lgdistplot(self.lgdist_rv_data, 'Reverse')
            )


    #Parse a generic substitution YAML file (3' and 5' supported)
    def parseFreqPlot(self, f, dict_to_add):

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
        
        #Sample name is always the only key in each of the supplied YAML files
        if parsed_data is not None:
            key = list(parsed_data.keys())[0] 
            s_name = self.clean_s_name(list(parsed_data.keys())[0],'')

            if s_name in dict_to_add:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            # Create tuples out of entries
            pos = list(range(1,len(parsed_data.get(key))))
            tuples = list(zip(pos,parsed_data.get(key)))
            # Get a dictionary out of it
            data = dict((x, y) for x, y in tuples)
            dict_to_add[s_name] = data
        else: 
            log.debug('No valid data {} in report'.format(f['fn']))
            return None

    #Parse a generic lgdistribution file and parse it to data frame
    def parselgDist(self, f, dict_to_add):

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
        
        #Sample name is always the only key in each of the supplied YAML files
        if parsed_data is not None:
            key = list(parsed_data.keys())[0] 
            s_name = self.clean_s_name(key,'')

            if s_name in dict_to_add:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            dict_to_add[s_name] = parsed_data.get(key)
        else: 
            log.debug('No valid data {} in report'.format(f['fn']))
            return None

    
    #### Tables from here on 
    def dmgprof_misinc_stats(self, dict_to_plot, title, description):
        """ Take the parsed stats from the DamageProfiler and add it to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['damageprofiler'] = {
            'title': title,
            'description': description,
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'RdOr',
            'format': '{:,.0f}',
            'modify': lambda x: x * 100.0
        }
        self.general_stats_addcols(dict_to_plot, headers)
    

    ##TODO add here Table info from lgdistribution 

    #### Plotting from here on
    #Nice Linegraph plot for lgdist data

    def lgdistplot(self,dict_to_use,orientation):
        """Generate a read length distribution plot"""

        data = dict()
        for s_name in dict_to_use:
            try:
                data[s_name] = {int(d): int (dict_to_use[s_name][d]) for d in dict_to_use[s_name]}
            except KeyError:
                pass
        if len(data) == 0:
            log.debug('No valid data for forward read lgdist input!')
            return None
        
        config = {
            'id': 'length-distribution-' + orientation,
            'smooth_points': 50,
            'title': 'DamageProfiler: Read length distribution: ' + orientation,
            'ylab': 'Number of reads',
            'xlab': 'Readlength (bp)',
            'tt_label': '{point.y} reads of length {point.x}',
            'ymin': 0,
            'xmin': 0
        }
        return linegraph.plot(data,config)


    #Linegraph plot for 3pGtoA
    def threeprime_plot(self):
        """Generate a 3' G -> A linegraph plot"""

        data = dict()
        for s_name in self.threepGtoAfreq_data:
            try:
                data[s_name] = {int(d): float (self.threepGtoAfreq_data[s_name][d])*100 for d in self.threepGtoAfreq_data[s_name]}
            except KeyError:
                pass
        if len(data) == 0:
            log.debug('No valid data for 3\' G -> A input!')
            return None

        config = {
            'id': 'threeprime_misinc_plot',
            'title': 'DamageProfiler: 3\' G -> A misincorporation plot',
            'ylab': '% G to A substituted',
            'xlab': 'Nucleotide position from 3\'',
            'tt_label': '{point.y:.2f} % G -> A misincorporations at nucleotide position {point.x}',
            'ymin': 0,
            'xmin': 1
        }

        return linegraph.plot(data,config)

    #Linegraph plot for 3pGtoA
    def fiveprime_plot(self):
        """Generate a 5' C -> T linegraph plot"""

        data = dict()
        for s_name in self.fivepCtoTfreq_data:
            try:
                data[s_name] = {int(d): float (self.fivepCtoTfreq_data[s_name][d])*100 for d in self.fivepCtoTfreq_data[s_name]}
            except KeyError:
                pass
        if len(data) == 0:
            log.debug('No valid data for 5\' C -> T input!')
            return None

        config = {
            'id': 'fiveprime_misinc_plot',
            'title': 'DamageProfiler: 5P C -> T misincorporation plot',
            'ylab': '% C to T substituted',
            'xlab': 'Nucleotide position from 5\'',
            'tt_label': '{point.y:.2f} % C -> T misincorporations at nucleotide position {point.x}',
            'ymin': 0,
            'xmin': 1
        }

        return linegraph.plot(data,config)

        
        