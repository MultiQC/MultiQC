""" MultiQC module to parse output from metadmg """


import logging
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .dfit import DfitReportMixin
from .stat import StatReportMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """ """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="metadmg",
            anchor="metadmg",
            href="https://github.com/metaDMG-dev/metaDMG-cpp",
            info="Taxonomic classification of reads and estimation of damage rates in ancient DNA data",
            doi="10.1101/2022.12.06.519264",
        )

        # Read config
        self.rank = getattr(config, "metadmg", {}).get("rank", "genus")
        self.n_taxa = getattr(config, "metadmg", {}).get("n_taxa", 10)
        self.index = getattr(config, "metadmg", {}).get("index", "name")
        self.sort_by = getattr(config, "metadmg", {}).get("sort_by", "nreads")

        # Set up class objects to hold parsed data
        n = dict()

        # Call submodule functions
        n["dfit"] = self.parse_metadmg_dfit()
        if n["dfit"] > 0:
            log.info("Found {} dfit reports".format(n["dfit"]))

        n["stat"] = self.parse_metadmg_stat()
        if n["stat"] > 0:
            log.info("Found {} stat reports".format(n["stat"]))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
