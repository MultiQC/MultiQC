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
from . import read_distribution
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
        self.sections = list()
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()
        
        # Call submodule functions
        n['bam_stat'] = bam_stat.parse_reports(self)
        if n['bam_stat'] > 0:
            log.info("Found {} bam_stat reports".format(n['bam_stat']))
        
        n['read_distribution'] = read_distribution.parse_reports(self)
        if n['read_distribution'] > 0:
            log.info("Found {} read_distribution reports".format(n['read_distribution']))
        
        n['gene_body_coverage'] = gene_body_coverage.parse_reports(self)
        if n['gene_body_coverage'] > 0:
            log.info("Found {} gene_body_coverage reports".format(n['gene_body_coverage']))
        
        n['inner_distance'] = inner_distance.parse_reports(self)
        if n['inner_distance'] > 0:
            log.info("Found {} inner_distance reports".format(n['inner_distance']))
        
        n['junction_annotation'] = junction_annotation.parse_reports(self)
        if n['junction_annotation'] > 0:
            log.info("Found {} junction_annotation reports".format(n['junction_annotation']))
        
        n['junction_saturation'] = junction_saturation.parse_reports(self)
        if n['junction_saturation'] > 0:
            log.info("Found {} junction_saturation reports".format(n['junction_saturation']))
        
        n['read_gc'] = read_gc.parse_reports(self)
        if n['read_gc'] > 0:
            log.info("Found {} read_gc reports".format(n['read_gc']))
        
        n['read_duplication'] = read_duplication.parse_reports(self)
        if n['read_duplication'] > 0:
            log.info("Found {} read_duplication reports".format(n['read_duplication']))
        

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning
        
        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)

    
