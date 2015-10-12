#!/usr/bin/env python

""" MultiQC module to parse output from Preseq """

from __future__ import print_function
import io
import logging
import os
import re

from multiqc import config, BaseMultiqcModule

# Initialise the logger
log = logging.getLogger()

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Preseq', anchor='preseq',
        href='http://smithlabresearch.org/software/preseq/', 
        info="estimates the complexity of a library, showing how many additional "\
         "unique reads are sequenced for increasing total read count.")

        # Find and load any Preseq reports
        self.preseq_data = dict()
        
        for f in self.find_log_files('ccurve.txt'):
            parsed_data = self.parse_preseq_logs(f)
            if parsed_data is not None:
                self.preseq_data[f['s_name']] = parsed_data

        if len(self.preseq_data) == 0:
            log.debug("Could not find any preseq data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.preseq_data)))

        self.sections = list()

        # Preseq plot
        # Only one section, so add to the intro
        self.intro += self.preseq_length_trimmed_plot()


    def parse_preseq_logs(self, f):
        """ Go through log file looking for preseq output """
        
        lines = f['f'].splitlines()
        header = lines.pop(0)
        if header != 'TOTAL_READS     EXPECTED_DISTINCT       LOWER_0.95CI    UPPER_0.95CI':
            log.debug("First line of preseq file {} did not look right".format(f['fn']))
            return None
        
        data = {}
        for l in lines:
            s = l.split()
            data[float(s[1])] = float(s[0])
        return data


    def preseq_length_trimmed_plot (self):
        """ Generate the preseq plot """    
        pconfig = {
            'id': 'preseq_plot',
            'title': 'Preseq complexity curve',
            'ylab': 'Unique Molecules',
            'xlab': 'Total Molecules (including duplicates)',
            'ymin': 0,
            'xmin': 0,
        }
        return self.plot_xy_data(self.preseq_data, pconfig)
