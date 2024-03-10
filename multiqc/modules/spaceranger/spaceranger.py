""" MultiQC module to parse output from Space Ranger """

import logging

from multiqc.modules.base_module import BaseMultiqcModule

# Import the Samtools submodules
from .count import SpaceRangerCountMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, SpaceRangerCountMixin):
    """Space Ranger has 2 main modules: count and vdj.
    This module parses data directly from the web_summary.html and summarise
    data useful for QC in the main sample table as well as generate some module specific plots
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Space Ranger",
            anchor="spaceranger",
            href="https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger",
            info="Space Ranger analyzes 10x Genomics Visium spatial transcriptomics data.",
            doi=[],
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = dict()
        self.general_stats_data = dict()
        n = dict()

        # Call submodule functions
        n["count"] = self.parse_count_html()
        if n["count"] > 0:
            log.info("Found {} Space Ranger count reports".format(n["count"]))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
