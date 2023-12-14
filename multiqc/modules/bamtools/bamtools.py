""" MultiQC module to parse output from Bamtools """

import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Import the Bamtools submodules
from . import stats

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """Bamtools is a collection of scripts. This MultiQC module
    supports some but not all. The code for each script is split
    into its own file and adds a section to the module ooutput if
    logs are found."""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Bamtools",
            anchor="bamtools",
            href="https://github.com/pezmaster31/bamtools",
            info="provides both a programmer's API and an end-user's toolkit for handling BAM files.",
            doi="10.1093/bioinformatics/btr174",
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = dict()
        self.general_stats_data = dict()
        n = dict()

        # Call submodule functions
        n["stats"] = stats.parse_reports(self)
        if n["stats"] > 0:
            log.info(f"Found {n['stats']} bamtools stats reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
