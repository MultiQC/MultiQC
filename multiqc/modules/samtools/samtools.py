""" MultiQC module to parse output from Samtools """


import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

from .stats import StatsReportMixin
from .flagstat import FlagstatReportMixin
from .idxstats import IdxstatsReportMixin
from .rmdup import RmdupReportMixin
from .coverage import CoverageReportMixin
from .markdup import MarkdupReportMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(
    BaseMultiqcModule,
    StatsReportMixin,
    FlagstatReportMixin,
    IdxstatsReportMixin,
    RmdupReportMixin,
    MarkdupReportMixin,
    CoverageReportMixin,
):
    """Samtools has a number of different commands and outputs.
    This MultiQC module supports some but not all. The code for
    each script is split into its own file and adds a section to
    the module output if logs are found."""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Samtools",
            anchor="Samtools",
            target="Samtools",
            href="http://www.htslib.org",
            info=" is a suite of programs for interacting with high-throughput sequencing data.",
            doi="10.1093/bioinformatics/btp352",
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = dict()
        self.general_stats_data = dict()
        n = dict()

        # Call submodule functions
        n["stats"] = self.parse_samtools_stats()
        if n["stats"] > 0:
            log.info(f"Found {n['stats']} stats reports")

        n["flagstat"] = self.parse_samtools_flagstat()
        if n["flagstat"] > 0:
            log.info(f"Found {n['flagstat']} flagstat reports")

        n["idxstats"] = self.parse_samtools_idxstats()
        if n["idxstats"] > 0:
            log.info(f"Found {n['idxstats']} idxstats reports")

        n["rmdup"] = self.parse_samtools_rmdup()
        if n["rmdup"] > 0:
            log.info(f"Found {n['rmdup']} rmdup reports")

        n["coverage"] = self.parse_samtools_coverage()
        if n["coverage"] > 0:
            log.info(f"Found {n['coverage']} coverage reports")

        n["markdup"] = self.parse_samtools_markdup()
        if n["markdup"] > 0:
            log.info(f"Found {n['markdup']} markdup reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
