""" MultiQC module to parse output from HOMER """


import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Import the HOMER submodules
from .findpeaks import FindPeaksReportMixin
from .tagdirectory import TagDirReportMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, FindPeaksReportMixin, TagDirReportMixin):
    """HOMER has a number of different commands and outputs.
    This MultiQC module supports some but not all. The code for
    each script is split into its own file and adds a section to
    the module output if logs are found."""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="HOMER",
            anchor="homer",
            href="http://homer.ucsd.edu/homer/",
            info="is a suite of tools for Motif Discovery and next-gen sequencing analysis.",
            doi="10.1016/j.molcel.2010.05.004",
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = dict()
        self.general_stats_data = dict()
        n = dict()

        # Set up data structures
        self.tagdir_data = {
            "GCcontent": {},
            "restriction": {},
            "restriction_norm": {},
            "length": {},
            "taginfo_total": {},
            "taginfo_total_norm": {},
            "taginfo_uniq": {},
            "taginfo_uniq_norm": {},
            "FreqDistribution": {},
            "header": {},
            "interChr": {},
        }
        # Call submodule functions

        n["findpeaks"] = self.parse_homer_findpeaks()
        if n["findpeaks"] > 0:
            log.info(f"Found {n['findpeaks']} findPeaks reports")

        n["tagDir"] = self.homer_tagdirectory()
        if n["tagDir"] > 0:
            log.info(f"Found {n['tagDir']} tagDir reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
