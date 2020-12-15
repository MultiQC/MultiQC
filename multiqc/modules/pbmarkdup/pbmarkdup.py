#!/usr/bin/env python

""" MultiQC module to parse output from pbmarkdup"""

import logging

from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ 
    pbmarkdup module class, parses pbmarkdup output.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
                name='pbmarkdup', anchor='pbmarkdup',
                href='https://github.com/PacificBiosciences/pbmarkdup',
                info=("""
                 takes one or multiple sequencing chips of an amplified libray
                 as HiFi reads and marks or removes duplicates.
                 """)
        )

        self.pbmarkdup = dict()

        for logfile in self.find_log_files('pbmarkdup', filehandles=True):
            print(logfile['fn'])
