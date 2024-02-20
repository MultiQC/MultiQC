""" MultiQC module to parse output from Lima """

import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Import HUMID submodules
from . import stats, neighbours, counts, clusters

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialse the parent object
        super(MultiqcModule, self).__init__(
            name="HUMID",
            anchor="humid",
            href="https://github.com/jfjlaros/HUMID",
            info=" is a fast, reference free tool to remove (UMI) duplicates from sequencing data",
            # No publication / DOI // doi=
        )
        self.stats = None
        self.neighbours = None
        self.counts = None
        self.clusters = None

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Look for stats files
        stats.parse_reports(self)

        # Look for neighbour files
        neighbours.parse_reports(self)

        # Look for count files
        counts.parse_reports(self)

        # Look for cluster files
        clusters.parse_reports(self)

        if all(not x for x in [self.stats, self.neighbours, self.counts, self.clusters]):
            raise ModuleNoSamplesFound

        num_samples = max([len(x) for x in [self.stats, self.neighbours, self.counts, self.clusters]])
        log.info(f"Found {num_samples} reports")
