""" MultiQC module to parse output from Sentieon-dnaseq """


import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule

# Import the Sentieon submodules
from . import AlignmentSummaryMetrics, GcBiasMetrics, InsertSizeMetrics

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """Sentieon-dnaseq produces many outputs. This module deals with 3 Picard
    equivalents which do not transfer well to MultiQC. The code for each script
    is split into its own file and adds a section to the module output if
    logs are found."""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Sentieon",
            anchor="sentieon",
            href="https://www.sentieon.com/products/",
            info="contains a suite of QC tools. The ones represented in this module are analogues of Picard QC metrics.",
            # Can't find a DOI // doi=
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

        # Call submodule functions
        n["AlignmentMetrics"] = AlignmentSummaryMetrics.parse_reports(self)
        if n["AlignmentMetrics"] > 0:
            log.info("Found {} AlignmentSummaryMetrics reports".format(n["AlignmentMetrics"]))

        n["GcBiasMetrics"] = GcBiasMetrics.parse_reports(self)
        if n["GcBiasMetrics"] > 0:
            log.info("Found {} GcBiasMetrics reports".format(n["GcBiasMetrics"]))

        n["InsertSizeMetrics"] = InsertSizeMetrics.parse_reports(self)
        if n["InsertSizeMetrics"] > 0:
            log.info("Found {} InsertSizeMetrics reports".format(n["InsertSizeMetrics"]))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise UserWarning

        # Add to the General Stats table (has to be called once per
        # MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)

    # Helper functions
    def multiply_hundred(self, val):
        try:
            val = float(val) * 100
        except ValueError:
            pass
        return val
