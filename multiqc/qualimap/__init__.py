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
             the quality control of alignment sequencing data and its derivatives like feature counts.</p>'
        self.analysis_dir = report['analysis_dir']
        self.output_dir = report['output_dir']
        self.data_dir = os.path.join(self.output_dir, 'report_data', 'qualimap')

        # Find QualiMap reports
        qualimap_raw_data = {}
        for root, dirnames, filenames in os.walk(self.analysis_dir, followlinks=True):
            if 'genome_results.txt' in filenames and 'raw_data_qualimapReport' in dirnames:
                with open(os.path.join(root, 'genome_results.txt'), 'r') as gr:
                    for l in gr:
                        if 'bam file' in l:
                            s_name = self.clean_s_name(os.path.basename(l.split(' = ')[-1]))

                if report['prepend_dirs']:
                    s_name = "{} | {}".format(os.path.dirname(root).replace(os.sep, ' | '), s_name).lstrip('. | ')

                qualimap_raw_data[s_name] = {'reports': [], 'plots': []}
                [qualimap_raw_data[s_name]['reports'].append(r) for r in \
                    os.listdir(os.path.join(root, 'raw_data_qualimapReport'))]
                [qualimap_raw_data[s_name]['plots'].append(p) for p in \
                    os.listdir(os.path.join(root, 'images_qualimapReport'))]

        if len(qualimap_raw_data) == 0:
            logging.debug("Could not find any QualiMap reports in {}".format(self.analysis_dir))
            raise UserWarning

        logging.info("Found {} QualiMap reports".format(len(qualimap_raw_data)))
