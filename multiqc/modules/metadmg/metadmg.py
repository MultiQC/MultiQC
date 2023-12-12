""" MultiQC module to plot output from metadmg """


import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .stat import StatReportMixin
from .dfit import DfitReportMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, StatReportMixin, DfitReportMixin):
    """
    To allow reading gzip archives, run with `ignore_images: false`
    in the config, e.g.:
    ```
    multiqc . --cl-config 'ignore_images: false'
    ```
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="metaDMG",
            anchor="metadmg",
            href="https://github.com/metaDMG-dev/metaDMG-cpp",
            info="Taxonomic classification of reads and estimation of damage rates in ancient DNA data",
            doi="10.1101/2022.12.06.519264",
        )

        # Read config
        self.rank = getattr(config, "metadmg", {}).get("rank", "genus")
        self.n_taxa = getattr(config, "metadmg", {}).get("n_taxa", 20)
        self.sort_by = getattr(config, "metadmg", {}).get("sort_by", "nreads")

        # Set up class objects to hold parsed data
        n = dict()

        # Call submodule functions
        n["stat"] = self.parse_metadmg_stat()
        if n["stat"] > 0:
            log.info("Found {} stat reports".format(n["stat"]))

        # n["dfit"] = self.parse_metadmg_dfit()
        # if n["dfit"] > 0:
        #    log.info("Found {} dfit reports".format(n["dfit"]))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
