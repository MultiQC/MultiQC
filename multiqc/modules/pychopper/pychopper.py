#!/usr/bin/env python

""" MultiQC module to parse output from pychopper """

from __future__ import print_function
import logging
import os
from collections import OrderedDict

from multiqc.plots import linegraph, bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Pychopper', anchor='pychopper',
                href='https://github.com/nanoporetech/pychopper',
        info="is a tool to identify, orient, trim and rescue full length Nanopore cDNA reads")
        
        # Parse stats file
        self.pychopper_data = {}
        for f in self.find_log_files('pychopper'):
            sample = f['s_name']
            self.pychopper_data[sample]= {}
            lines = f['f'].splitlines()
            for line in lines[1:]:
                category, name, value = line.split()
                if category not in self.pychopper_data[sample]:
                    self.pychopper_data[sample][category] = {}
                self.pychopper_data[sample][category][name] = float(value)
         
        # Add to general statistics table:
        # Percentage of full length transcripts
        data_general_stats={}
        for sample in self.pychopper_data.keys():
            data_general_stats[sample] = {}
            c=self.pychopper_data[sample]['Classification']
            ftp=c['Primers_found'] * 100 / (c['Primers_found'] + c['Rescue'] + c['Unusable'])
            data_general_stats[sample]['ftp']=ftp

        headers = OrderedDict()
        headers['ftp'] = {
            'title' : 'Full length Transp.',
            'description' : 'Percentage of full length transposons in the cDNA reads',
            'suffix' : '%',
            'max': 100,
            'min': 0,
        }
        
        self.general_stats_addcols(data_general_stats, headers)

        # Write data file
        self.write_data_file(self.pychopper_data, 'multiqc_pychopper')

        # Report sections
        self.add_section (
            name = "cDNA Read Classification",
            description = (
                """
                This plot shows the read categories identified by porechopper: 
                    - Full length transcripts (primers found at both ends),
                    - Rescued reads (Splitted fusion reads)
                    - Unusable (No primer sequences on both ends)
                """
            ),
            helptext = (
                """
                TODO: Add helptext here
                """
            ),
            anchor = 'pychopper_classification',
            plot = self.plot_classification()
        )
        self.add_section (
            name = "cDNA Strand Orientation",
            description = (
                """
                This plot shows the idientfied orientation of full length transcripts
                """
            ),
            helptext = (
                """
                TODO: Add helptext here
                """
            ),
            anchor = 'pychopper_orientation',
            plot = self.plot_orientation()
        )

    # Plotting functions
    def plot_classification(self):
        """ Generate the cDNA read classification plot """

        pconfig = {
            'id': 'pychopper_classification',
            'title': 'Pychopper: Read classification'
     #       'ylab': 'Sample',
     #       'xlab': 'X-Axis',
     #       'xDecimals': False,
     #       'ymin': 0
        }

        data_classification = {}
        for sample in self.pychopper_data.keys():
            data_classification[sample] = {}
            data_classification[sample]=self.pychopper_data[sample]['Classification']
        
        print(data_classification)
        cats = ['Primers_found', 'Rescue', 'Unusable']
        return bargraph.plot(data_classification, cats, pconfig)
    
    def plot_orientation(self):
        """ Generate the transcript strand orientation plot """

        pconfig = {
            'id': 'pychopper_orientation',
            'title': 'Pychopper: Strand Orientation',
            'ylab': 'Count',
            'xlab': 'Values',
            'xDecimals': False,
            'ymin': 0
        }

        data_orientation = {}
        for sample in self.pychopper_data.keys():
            data_orientation[sample] = {}
            data_orientation[sample]=self.pychopper_data[sample]['Strand']

        print(data_orientation)
        cats = ['+', '-']
        return bargraph.plot(data_orientation, cats, pconfig)

        

