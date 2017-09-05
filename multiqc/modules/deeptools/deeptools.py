#!/usr/bin/env python
"""MultiQC module to parse the output from deepTools"""
from collections import OrderedDict
import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

# deepTools modules
from .bamPEFragmentSize import bamPEFragmentSizeMixin
from .estimateReadFiltering import estimateReadFilteringMixin
from .plotCoverage import plotCoverageMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, bamPEFragmentSizeMixin, estimateReadFilteringMixin, plotCoverageMixin):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='deepTools', anchor='deepTools', target='deepTools',
                                            href="http://deeptools.readthedocs.io",
                                            info=" is a suite of tools to process and analyze deep sequencing data.")

        # Set up class objects to hold parsed data
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

        # bamPEFragmentSize
        n['bamPEFragmentSize'] = self.parse_bamPEFragmentSize()
        if n['bamPEFragmentSize'] > 0:
            log.info("Found {} deepTools 'bamPEFragmentSize --table' samples".format(n['bamPEFragmentSize']))
        else:
            log.debug("Could not find any 'bamPEFragmentSize --table' outputs in {}".format(config.analysis_dir))

        # estimateReadFiltering
        n['estimateReadFiltering'] = self.parse_estimateReadFiltering()
        if n['estimateReadFiltering'] > 0:
            log.info("Found {} deepTools estimateReadFiltering samples".format(n['estimateReadFiltering']))
        else:
            log.debug("Could not find any estmateReadFiltering outputs in {}".format(config.analysis_dir))

        # plotCoverage
        n['plotCoverageStdout'], n['plotCoverageOutRawCounts'] = self.parse_plotCoverage()
        if n['plotCoverageStdout'] + n['plotCoverageOutRawCounts'] > 0:
            extra = ''
            if n['plotCoverageOutRawCounts'] == 0:
                extra = ' (you may need to increase the maximum log file size to find plotCoverage --outRawCounts files)'
            log.info("Found {} and {} deepTools plotCoverage standard output and --outRawCounts samples, respectively{}".format(n['plotCoverageStdout'], n['plotCoverageOutRawCounts'], extra))
        else:
            log.debug("Could not find any plotCoverage outputs in {} (you may need to increase the maximum log file size)".format(config.analysis_dir))
