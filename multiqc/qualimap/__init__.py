#!/usr/bin/env python

""" MultiQC module to parse output from QualiMap """

from __future__ import print_function
import logging
import os

import multiqc

class MultiqcModule(multiqc.BaseMultiqcModule):

    def __init__(self, report):

        # Initialise the parent object
        super(MultiqcModule, self).__init__()

        # Static variables
        self.name = "QualiMap"
        self.anchor = "qualimap"
        self.intro = '<p><a href="http://qualimap.bioinfo.cipf.es/" target="_blank">QualiMap</a> \
             is a platform-independent application written in Java and R that provides both \
			 a Graphical User Inteface (GUI) and a command-line interface to facilitate \
			 the quality control of alignment sequencing data and its derivatives like feature counts. .</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']
