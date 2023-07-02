"""MultiQC module to parse the output from sincei"""
import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule

# sincei modules
from .scCountQC import scCountQCMixin
from .scFilterStats import scFilterStatsMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(
    BaseMultiqcModule,
    scFilterStatsMixin,
    scCountQCMixin,
):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="sincei",
            anchor="sincei",
            target="sincei",
            href="http://sincei.readthedocs.io",
            info=" is a toolkit to process and analyze single-cell (epi)genomics data.",
            doi="10.5281/zenodo.7853375",
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

        # scFilterStats
        n["scFilterStats"] = self.parse_scFilterStats()
        if n["scFilterStats"] > 0:
            log.debug("Found {} sincei scFilterStats samples".format(n["scFilterStats"]))

        # scCountQC
        n["scCountQC"] = self.parse_scCountQC()
        if n["scCountQC"] > 0:
            log.debug("Found {} sincei scCountQC samples".format(n["scCountQC"]))

        tot = sum(n.values())
        if tot > 0:
            log.info("Found {} total sincei samples".format(tot))
        else:
            raise UserWarning

    def _int(self, val):
        """Avoids Python3 error:
        ValueError: invalid literal for self._int() with base 10: '1.0'
        """
        return int(round(float(val)))
