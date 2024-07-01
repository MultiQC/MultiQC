"""MultiQC module to parse output from bcftools"""

import logging
from typing import Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Import the Samtools submodules
from multiqc.modules.bcftools.stats import parse_bcftools_stats

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
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
        self.general_stats_headers: Dict = dict()
        self.general_stats_data: Dict = dict()
        n = dict()

        # Call submodule functions
        n["stats"] = parse_bcftools_stats(self)
        if n["stats"] > 0:
            log.info(f"Found {n['stats']} stats reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
