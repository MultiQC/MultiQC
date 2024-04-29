"""MultiQC submodule to parse output from Glimpse concordance analysis"""

import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound


from .err_spl import ErrSplReportMixin  # Import the Glimpse submodules

# Initialise the logger
log = logging.getLogger(__name__)


# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, ErrSplReportMixin):
    """Glimpse has a number of different commands and outputs.
    This MultiQC module supports some but not all."""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Glimpse",
            anchor="glimpse",
            target="Glimpse",
            href="https://odelaneau.github.io/GLIMPSE/",
            info="Set of tools for low-coverage whole genome sequencing imputation",
            doi="10.1101/2022.11.28.518213 ",
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = dict()
        self.general_stats_data = dict()
        n = dict()

        # Call submodule functions
        n["err_spl"] = self.parse_err_spl()
        if n["err_spl"] > 0:
            log.info(f"Found {n['stats']} errors by sample reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
