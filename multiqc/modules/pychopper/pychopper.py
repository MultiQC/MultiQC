#!/usr/bin/env python

""" MultiQC module to parse output from pychopper """

from __future__ import print_function
from collections import OrderedDict
import logging
import os

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

        print("This is the pychopper Module writing")
        
        self.pychopper_data = dict()
        for f in self.find_log_files('pychopper'):
            print(f)
        
        # Add to general statistics table


        # Write data file
        self.write_data_file(self.pychopper_data, 'multiqc_pychopper')

        # Report sections
        self.add_section (
            name = "cDNA Read Classification"
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
            anchor = 'pychopper-classification',
            plot = self.plot_classification()
        )
        self.add_section (
            name = "cDNA Strand Orientation"
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
            anchor = 'pychopper-orientation',
            plot = self.plot_orientation()
        )

    # Plotting functions
    def plot_classification(self):
        """ Generate the cDNA read classification plot """

    pconfig = {
        'id': 'pychopper_classification',
        'title': 'Pychopper: Read classification',
        'ylab': 'Sample',
        'xlab': 'X-Axis',
        'xDecimals': False,
        'ymin': 0
    }

    return bargraph.plot(self.pychopper_classification, pconfig)
    

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

    return bargraph.plot(self.pychopper_orientation, pconfig)

