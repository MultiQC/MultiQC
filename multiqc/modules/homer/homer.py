#!/usr/bin/env python
""" MultiQC module to parse output from HOMER """
from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

# Import the HOMER submodules
from .findpeaks import FindPeaksReportMixin
from .tagdirectory import TagDirReportMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, FindPeaksReportMixin, TagDirReportMixin):
    """ HOMER has a number of different commands and outputs.
    This MultiQC module supports some but not all. The code for
    each script is split into its own file and adds a section to
    the module output if logs are found. """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='HOMER', anchor='homer',
            href='http://homer.ucsd.edu/homer/',
            info="is a suite of tools for Motif Discovery and next-gen sequencing analysis.")


        # Set up class objects to hold parsed data
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

         # Set up data structures
        self.tagdir_data = {
            'GCcontent': {},
            'restriction': {},
            'restriction_norm': {},
            'length': {},
            'taginfo_total': {},
            'taginfo_total_norm': {},
            'taginfo_uniq': {},
            'taginfo_uniq_norm': {},
            'FreqDistribution' : {},
            'header' : {},
            'interChr': {}

        }
        # Call submodule functions
        

        n['findpeaks'] = self.parse_homer_findpeaks()
        if n['findpeaks'] > 0:
            log.info("Found {} findPeaks reports".format(n['findpeaks']))

    
        n['Homer_tagDir'] = self.homer_tagdirectory()
        if n['Homer_tagDir'] > 0:
            log.info("Found {} reports".format(n['Homer_tagDir']))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
