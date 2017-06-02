#!/usr/bin/env python

""" MultiQC module to parse output from RSeQC """

from collections import OrderedDict
import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

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
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

        # Get the list of submodules (can be customised)
        rseqc_sections = getattr(config, 'rseqc_sections', [])
        if len(rseqc_sections) == 0:
            rseqc_sections = [
                'read_distribution',
                'gene_body_coverage',
                'inner_distance',
                'read_gc',
                'read_duplication',
                'junction_annotation',
                'junction_saturation',
                'infer_experiment',
                'bam_stat'
            ]

        # Call submodule functions
        for sm in rseqc_sections:
            try:
                # Import the submodule and call parse_reports()
                #   Function returns number of parsed logs
                module = __import__('multiqc.modules.rseqc.{}'.format(sm), fromlist=[''])
                n[sm] = getattr(module, 'parse_reports')(self)
                if n[sm] > 0:
                    log.info("Found {} {} reports".format(n[sm], sm))
            except (ImportError, AttributeError):
                log.warn("Could not find RSeQC Section '{}'".format(sm))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)


