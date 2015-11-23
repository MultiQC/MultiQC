#!/usr/bin/env python

""" MultiQC module to parse output from QualiMap """

from __future__ import print_function
from collections import OrderedDict
import io
import logging
import os

from collections import defaultdict

from multiqc import config, BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ Qualimap is really a collection of separate programs:
    BamQC, RNASeq and Counts.. This module is split into separate
    files to reflect this and help with code organisation. """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='QualiMap', anchor='qualimap',
        href="http://qualimap.bioinfo.cipf.es/", 
        info="is a platform-independent application to facilitate the quality"\
        " control of alignment sequencing data and its derivatives like"\
        " feature counts.")
        
        # Global dict used by all submodules
        self.general_stats = defaultdict(dict)
        
        # Initialise the BamQC submodule and parse logs
        import BamQC
        BamQC.parse_reports(self)


        if len(self.general_stats) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.general_stats)))

        # Add the submodule results to the report output
        self.sections = list()
        BamQC.report_sections(self)

        # General stats table columns
        BamQC.stats_table(self)

    
