""" MultiQC module to parse output from Sambamba """


import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Import the Sambamba submodules
from .markdup import SambambaMarkdupMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, SambambaMarkdupMixin):
    """Sambamba has a number of different commands and outputs.
    This MultiQC module supports some but not all. The code for
    each script is split into its own file and adds a section to
    the module output if logs are found."""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Sambamba",
            anchor="sambamba",
            href="https://lomereiter.github.io/sambamba/",
            info=" is a suite of programs for interacting with high-throughput sequencing data.",
            doi="10.1093/bioinformatics/btv098",
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = dict()
        self.general_stats_data = dict()
        n = dict()

        # Call submodule functions
        n["markdup"] = self.parse_sambamba_markdup()

        if n["markdup"] > 0:
            log.info(f"Found {n['markdup']} sambamba markdup reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
