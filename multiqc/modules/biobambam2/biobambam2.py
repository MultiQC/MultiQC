""" MultiQC module to parse output from biobambam2 """


import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.modules.picard import MarkDuplicates

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """This module is super weird. The output from this tools is essentially
    identical to Picard MarkDuplicates, so we just hijack that module instead"""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="biobambam2",
            anchor="biobambam2",
            href="https://gitlab.com/german.tischler/biobambam2",
            info="provides tools for early stage alignment file processing",
            doi="10.1186/1751-0473-9-13",
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = dict()
        self.general_stats_data = dict()
        n = dict()

        n["bamsormadup"] = MarkDuplicates.parse_reports(self, "biobambam2/bamsormadup")
        if n["bamsormadup"] > 0:
            log.info(f"Found {n['bamsormadup']} bamsormadup reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)

    # Helper functions
    @staticmethod
    def multiply_hundred(val):
        try:
            val = float(val) * 100
        except ValueError:
            pass
        return val
