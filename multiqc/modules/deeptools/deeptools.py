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
from .plotEnrichment import plotEnrichmentMixin
from .plotFingerprint import plotFingerprintMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, bamPEFragmentSizeMixin, estimateReadFilteringMixin, plotCoverageMixin, plotEnrichmentMixin, plotFingerprintMixin):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='deepTools', anchor='deepTools', target='deepTools',
                                            href="http://deeptools.readthedocs.io",
                                            info=" is a suite of tools to process and analyze deep sequencing data.")

        # Set up class objects to hold parsed data
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

        # estimateReadFiltering
        n['estimateReadFiltering'] = self.parse_estimateReadFiltering()
        if n['estimateReadFiltering'] > 0:
            log.debug("Found {} deepTools estimateReadFiltering samples".format(n['estimateReadFiltering']))

        # plotCoverage
        n['plotCoverageStdout'], n['plotCoverageOutRawCounts'] = self.parse_plotCoverage()
        if n['plotCoverageStdout'] + n['plotCoverageOutRawCounts'] > 0:
            extra = ''
            if n['plotCoverageOutRawCounts'] == 0:
                extra = ' (you may need to increase the maximum log file size to find plotCoverage --outRawCounts files)'
            log.debug("Found {} and {} deepTools plotCoverage standard output and --outRawCounts samples, respectively{}".format(n['plotCoverageStdout'], n['plotCoverageOutRawCounts'], extra))

        # bamPEFragmentSize
        n['bamPEFragmentSize'] = self.parse_bamPEFragmentSize()
        if n['bamPEFragmentSize'] > 0:
            log.debug("Found {} deepTools 'bamPEFragmentSize --table' samples".format(n['bamPEFragmentSize']))

        # plotEnrichment
        n['plotEnrichment'] = self.parse_plotEnrichment()
        if n['plotEnrichment'] > 0:
            log.debug("Found {} deepTools plotEnrichment samples".format(n['plotEnrichment']))

        # plotFingerprint
        n['plotFingerprintOutQualityMetrics'], n['plotFingerprintOutRawCounts'] = self.parse_plotFingerprint()
        if n['plotFingerprintOutQualityMetrics'] + n['plotFingerprintOutRawCounts'] > 0:
            extra = ''
            if n['plotFingerprintOutRawCounts'] == 0:
                extra = ' (you may need to increase the maximum log file size to find plotFingerprint --outRawCounts files)'
            log.debug("Found {} and {} deepTools plotFingerprint --outQualityMetrics and --outRawCounts samples, respectively{}".format(n['plotFingerprintOutQualityMetrics'], n['plotFingerprintOutRawCounts'], extra))

        tot = sum(n.values())
        if tot > 0:
            log.info('Found {} total deepTools samples'.format(tot))
        else:
            raise UserWarning
