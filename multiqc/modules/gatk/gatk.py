#!/usr/bin/env python
""" MultiQC module to parse output from GATK """
from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config, BaseMultiqcModule

# Import the GATK submodules
# import varianteval
from .varianteval import VariantEvalMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, VariantEvalMixin):
    """ GATK has a number of different commands and outputs.
    This MultiQC module supports some but not all. The code for
    each script is split into its own file and adds a section to
    the module output if logs are found. """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='GATK',  anchor='gatk', target='GATK',
            href='https://www.broadinstitute.org/gatk/',
            info=(" is a toolkit offering a wide variety of tools with a "
                  "primary focus on variant discovery and genotyping."))

        # Set up class objects to hold parsed data
        self.sections = list()
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

        # Call submodule functions
        n['varianteval'] = self.parse_gatk_varianteval()
        if n['varianteval'] > 0:
            log.info("Found {} VariantEval reports".format(n['varianteval']))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
