""" MultiQC submodule to parse output from Sentieon AlignmentSummaryMetrics
 (based on the Picard module of the same name) """

import logging

from multiqc.modules.picard import AlignmentSummaryMetrics

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """
    Sentieon uses Picard tools to generate metrics, so reusing the function
    from Picard tools to parse the reports.
    """

    return AlignmentSummaryMetrics.parse_reports(self)
