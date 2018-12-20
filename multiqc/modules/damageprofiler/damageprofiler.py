#!/usr/bin/env python

""" MultiQC module to parse output from DamageProfiler """

from __future__ import print_function
from collections import OrderedDict
import logging
import json

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

        # Init empty dictionaries

        self.threepGtoAfreq_data = dict()
        self.fivepCtoTfreq_data = dict()
        self.lgdist_fw_data = dict()
        self.lgdist_rv_data = dict()

        # Find and load JSON file
        for f in self.find_log_files('damageprofiler',filehandles=True):
            self.parseJSON(f)

        # Filter to strip out ignored sample names
        self.threepGtoAfreq_data         =   self.ignore_samples(self.threepGtoAfreq_data)
        self.fivepCtoTfreq_data          =   self.ignore_samples(self.fivepCtoTfreq_data)
        self.lgdist_fw_data   =   self.ignore_samples(self.lgdist_fw_data)
        self.lgdist_rv_data      =   self.ignore_samples(self.lgdist_rv_data)


        # Write parsed data to a file
        #self.write_data_file(self.threepGtoAfreq_data, 'multiqc_damageprofiler_3pGtoAfreq')
        #self.write_data_file(self.fivepCtoTfreq_data, 'multiqc_damageprofiler_5pCtoTfreq')
        #self.write_data_file(self.lgdist_fw_data, 'multiqc_damageprofiler_lgdist_fw')
        #self.write_data_file(self.lgdist_rv_data, 'multiqc_damageprofiler_lgdist_rv')

        # Basic Stats Table, use generic function to add data to general table
        #self.dmgprof_misinc_stats(self.threepGtoAfreq_data, '3 Prime', 'G -> A')
        #self.dmgprof_misinc_stats(self.fivepCtoTfreq_data, '5 Prime', 'C -> T')

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
              #  plot = self.lgdistplot(self.lgdist_fw_data, 'Forward')
            )
            self.add_section (
                name = 'Reverse read length distribution',
              #  plot = self.lgdistplot(self.lgdist_rv_data, 'Reverse')
            )


    #Parse our nice little JSON file
    def parseJSON(self, f):
        """ Parse the JSON output from DamageProfiler and save the summary statistics """
        try:
            parsed_json = json.load(f['f'])
        except Exception as e:
            print(e)
            log.warn("Could not parse DamageProfiler JSON: '{}'".format(f['fn']))
            return None
        
        #Get sample name from JSON first
        s_name = parsed_json['metadata']['sample_name']
        self.add_data_source(f, s_name)

        
        #Add 3' G to A data 
        self.threepGtoAfreq_data[s_name] = parsed_json['dmg_3p']

        #Add 5' C to T data
        self.fivepCtoTfreq_data[s_name] = parsed_json['dmg_5p']

        #Add lendist forward 
        self.lgdist_fw_data[s_name] = parsed_json['lendist_fw']

        #Add lendist reverse
        self.lgdist_rv_data[s_name] = parsed_json['lendist_rv']

    
    #### Tables from here on 

    def dmgprof_misinc_stats(self, dict_to_plot, readend, substitution):
        """ Take the parsed stats from the DamageProfiler and add it to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['1'] = {
            'title': readend + ' ' + substitution + ' 1st base ',
            'description': readend + ' ' + substitution + ' 1st base ',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'GrRd',
            'format': '{:,.2f}',
            'modify': lambda x: x * 100.0
        }
        headers['2'] = {
            'title': readend + ' ' + substitution + ' 2nd base ',
            'description': readend + ' ' + substitution + ' 2nd base ',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'GrRd',
            'format': '{:,.2f}',
            'modify': lambda x: x * 100.0
        }

        # Create new small subset dictionary for entries (we need just the first two data (k,v) pairs from each report)
        data = OrderedDict()
        
        for key in dict_to_plot.keys():
            tmp = dict_to_plot[key]
            #Extract first two elements from list
            data[key] = tmp[:2]
        
        self.general_stats_addcols(data,headers)
    

     

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
        dict_to_add = dict()
        # Create tuples out of entries
        for key in self.threepGtoAfreq_data: 
            pos = list(range(1,len(self.threepGtoAfreq_data.get(key))))
            tuples = list(zip(pos,self.threepGtoAfreq_data.get(key)))
            # Get a dictionary out of it
            data = dict((x, y) for x, y in tuples)
            dict_to_add[key] = data

        config = {
            'id': 'threeprime_misinc_plot',
            'title': 'DamageProfiler: 3\' G -> A misincorporation plot',
            'ylab': '% G to A substituted',
            'xlab': 'Nucleotide position from 3\'',
            'tt_label': '{point.y:.2f} % G -> A misincorporations at nucleotide position {point.x}',
            'ymin': 0,
            'xmin': 1
        }

        return linegraph.plot(dict_to_add,config)

    #Linegraph plot for 5pCtoT
    def fiveprime_plot(self):
        """Generate a 5' C -> T linegraph plot"""
        
        data = dict()
        dict_to_add = dict()
        # Create tuples out of entries
        for key in self.fivepCtoTfreq_data: 
            pos = list(range(1,len(self.fivepCtoTfreq_data.get(key))))
            tuples = list(zip(pos,self.fivepCtoTfreq_data.get(key)))
            # Get a dictionary out of it
            data = dict((x, y) for x, y in tuples)
            dict_to_add[key] = data

        config = {
            'id': 'fiveprime_misinc_plot',
            'title': 'DamageProfiler: 5\' C -> T misincorporation plot',
            'ylab': '% C to T substituted',
            'xlab': 'Nucleotide position from 5\'',
            'tt_label': '{point.y:.2f} % C -> T misincorporations at nucleotide position {point.x}',
            'ymin': 0,
            'xmin': 1
        }

        return linegraph.plot(dict_to_add,config)

        
        