#!/usr/bin/env python

""" MultiQC module to parse output from Picard """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config, BaseMultiqcModule

# Import the Picard submodules
from . import MarkDuplicates
from . import InsertSizeMetrics
from . import GcBiasMetrics
from . import HsMetrics
from . import OxoGMetrics

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ Picard is a collection of scripts. This MultiQC module
    supports some but not all. The code for each script is split
    into its own file and adds a section to the module output if
    logs are found."""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Picard', anchor='picard', 
        href='http://broadinstitute.github.io/picard/', 
        info="is a set of Java command line tools for manipulating high-"\
        "throughput sequencing data.")
        
        # Set up class objects to hold parsed data
        self.sections = list()
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()
        
        # Call submodule functions
        n['MarkDuplicates'] = MarkDuplicates.parse_reports(self)
        if n['MarkDuplicates'] > 0:
            log.info("Found {} MarkDuplicates reports".format(n['MarkDuplicates']))
        
        n['InsertSizeMetrics'] = InsertSizeMetrics.parse_reports(self)
        if n['InsertSizeMetrics'] > 0:
            log.info("Found {} InsertSizeMetrics reports".format(n['InsertSizeMetrics']))
        
        n['GcBiasMetrics'] = GcBiasMetrics.parse_reports(self)
        if n['GcBiasMetrics'] > 0:
            log.info("Found {} GcBiasMetrics reports".format(n['GcBiasMetrics']))
        
        n['HsMetrics'] = HsMetrics.parse_reports(self)
        if n['HsMetrics'] > 0:
            log.info("Found {} HsMetrics reports".format(n['HsMetrics']))
        
        n['OxoGMetrics'] = OxoGMetrics.parse_reports(self)
        if n['OxoGMetrics'] > 0:
            log.info("Found {} OxoGMetrics reports".format(n['OxoGMetrics']))
        
        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning
        
        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
    
    # Helper functions
    def multiply_hundred(self, val):
        try:
            val = float(val) * 100
        except ValueError:
            pass
        return val

