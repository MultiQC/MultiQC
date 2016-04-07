#!/usr/bin/env python

""" MultiQC module to parse output from RSeQC """

from collections import OrderedDict
import logging

from multiqc import config, BaseMultiqcModule

# Import the RSeQC submodules
from . import bam_stat
from . import gene_body_coverage
from . import inner_distance
from . import junction_annotation
from . import junction_saturation
from . import read_gc
# from . import read_distribution
from . import read_duplication

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ RSeQC is a collection of scripts. This MultiQC module
    supports some but not all. The code for each script is split
    into its own file and adds a section to the module ooutput if
    logs are found."""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='RSeQC', anchor='rseqc',
        href="http://rseqc.sourceforge.net/", 
        info="package provides a number of useful modules that can"\
        " comprehensively evaluate high throughput RNA-seq data.")
        
        # Set up class objects to hold parsed data
        self.sample_count = 0
        self.sections = list()
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        
        # Call submodule functions
        bam_stat.parse_reports(self)
        gene_body_coverage.parse_reports(self)
        inner_distance.parse_reports(self)
        junction_annotation.parse_reports(self)
        junction_saturation.parse_reports(self)
        read_gc.parse_reports(self)
        # read_distribution.parse_reports(self)
        read_duplication.parse_reports(self)

        # Exit if we didn't find anything
        if self.sample_count == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning
        
        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)

    
