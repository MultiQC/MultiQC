""" MultiQC submodule to parse output from Sentieon GcBiasMetrics
 (based on the Picard module of the same name """

import logging

from multiqc.modules.picard import GcBiasMetrics
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """
    Sentieon uses Picard tools to generate metrics, so reusing the function
    from Picard tools to parse the reports.
    """

    return GcBiasMetrics.parse_reports(self)
