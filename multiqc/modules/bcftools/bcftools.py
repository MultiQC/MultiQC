""" MultiQC module to parse output from bcftools """


import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Import the Samtools submodules
from .stats import StatsReportMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, StatsReportMixin):
    """Bcftools has a number of different commands and outputs.
    This MultiQC module supports some but not all. The code for
    each script is split into its own file and adds a section to
    the module output if logs are found."""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Bcftools",
            anchor="bcftools",
            target="Bcftools",
            href="https://samtools.github.io/bcftools/",
            info=" contains utilities for variant calling and manipulating VCFs and BCFs.",
            doi="10.1093/gigascience/giab008",
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

        # Call submodule functions
        n["stats"] = self.parse_bcftools_stats()
        if n["stats"] > 0:
            log.info("Found {} stats reports".format(n["stats"]))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
