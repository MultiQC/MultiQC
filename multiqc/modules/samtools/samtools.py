#!/usr/bin/env python

""" MultiQC module to parse output from Samtools """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config, BaseMultiqcModule

# Import the Samtools submodules
from . import stats

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ Samtools has a number of different commands and outputs.
    This MultiQC module supports some but not all. The code for
    each script is split into its own file and adds a section to
    the module output if logs are found. """
    
    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Samtools',
        anchor='Samtools', target='Samtools',
        href='http://www.htslib.org',
        info=" is a suite of programs for interacting with high-throughput sequencing data.")
        
        # Set up class objects to hold parsed data
        self.sections = list()
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()
        
        # Call submodule functions
        n['stats'] = stats.parse_reports(self)
        if n['stats'] > 0:
            log.info("Found {} stats reports".format(n['stats']))
        
        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning
        
        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
