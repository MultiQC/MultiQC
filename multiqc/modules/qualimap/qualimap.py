#!/usr/bin/env python

""" MultiQC module to parse output from QualiMap """

from __future__ import print_function
from collections import defaultdict, OrderedDict
import logging
import os

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

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

        # Initialise the submodules
        from . import QM_BamQC
        from . import QM_RNASeq

        # Set up class objects to hold parsed data()
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = defaultdict(lambda: dict())
        n = dict()

        # Call submodule functions
        n['BamQC'] = QM_BamQC.parse_reports(self)
        if n['BamQC'] > 0:
            log.info("Found {} BamQC reports".format(n['BamQC']))

        n['RNASeq'] = QM_RNASeq.parse_reports(self)
        if n['RNASeq'] > 0:
            log.info("Found {} RNASeq reports".format(n['RNASeq']))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)

    # Helper functions
    def get_s_name(self, f):
        s_name = os.path.basename(os.path.dirname(f['root']))
        s_name = self.clean_s_name(s_name, f['root'])
        if s_name.endswith('.qc'):
            s_name = s_name[:-3]
        return s_name
