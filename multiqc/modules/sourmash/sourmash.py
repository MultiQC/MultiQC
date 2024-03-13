""" MultiQC module to parse output from sourmash """

import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .compare import CompareMixin
from .gather import GatherMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, CompareMixin, GatherMixin):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Sourmash",
            anchor="sourmash",
            href="https://github.com/sourmash-bio/sourmash",
            info="quickly searches, compares, and analyzes genomic and metagenomic data sets.",
            doi="10.21105/joss.00027",
        )

        n = dict()
        n["compare"] = self.parse_compare()
        if n["compare"] > 0:
            log.info(f"Found {n['compare']} compare results")

        n["gather"] = self.parse_gather()
        if n["gather"] > 0:
            log.info(f"Found {n['gather']} gather results")

        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound
