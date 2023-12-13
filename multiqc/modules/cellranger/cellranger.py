""" MultiQC module to parse output from Cell Ranger """

import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Import the Samtools submodules
from .count import CellRangerCountMixin
from .vdj import CellRangerVdjMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, CellRangerCountMixin, CellRangerVdjMixin):
    """Cell Ranger has 2 main modules: count and vdj.
    This module parses data directly from the web_summary.html and summarise
    data useful for QC in the main sample table as well as generate some module specific plots
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Cell Ranger",
            anchor="cellranger",
            href="https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger",
            info="Cell Ranger analyzes single cell expression or VDJ data produced by 10X Genomics.",
            doi="10.1038/ncomms14049",
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = dict()
        self.general_stats_data = dict()
        n = dict()

        # Call submodule functions
        n["count"] = self.parse_count_html()
        if n["count"] > 0:
            log.info(f"Found {n['count']} Cell Ranger count reports")

        n["vdj"] = self.parse_vdj_html()
        if n["vdj"] > 0:
            log.info(f"Found {n['vdj']} Cell Ranger VDJ reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
