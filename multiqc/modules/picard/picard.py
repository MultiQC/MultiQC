"""MultiQC module to parse output from Picard"""

import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Import the Picard submodules, each one matching a picard tool
from . import (
    AlignmentSummaryMetrics,
    BaseDistributionByCycleMetrics,
    IlluminaBasecallingMetrics,
    IlluminaLaneMetrics,
    CrosscheckFingerprints,
    ExtractIlluminaBarcodes,
    GcBiasMetrics,
    HsMetrics,
    InsertSizeMetrics,
    MarkDuplicates,
    MarkIlluminaAdapters,
    OxoGMetrics,
    QualityByCycleMetrics,
    QualityScoreDistributionMetrics,
    QualityYieldMetrics,
    RnaSeqMetrics,
    RrbsSummaryMetrics,
    TargetedPcrMetrics,
    ValidateSamFile,
    VariantCallingMetrics,
    WgsMetrics,
)

TOOLS = [
    m.__name__.split(".")[-1]
    for m in (
        AlignmentSummaryMetrics,
        BaseDistributionByCycleMetrics,
        IlluminaBasecallingMetrics,
        IlluminaLaneMetrics,
        CrosscheckFingerprints,
        ExtractIlluminaBarcodes,
        GcBiasMetrics,
        HsMetrics,
        InsertSizeMetrics,
        MarkDuplicates,
        MarkIlluminaAdapters,
        OxoGMetrics,
        QualityByCycleMetrics,
        QualityScoreDistributionMetrics,
        QualityYieldMetrics,
        RnaSeqMetrics,
        RrbsSummaryMetrics,
        TargetedPcrMetrics,
        ValidateSamFile,
        VariantCallingMetrics,
        WgsMetrics,
    )
]

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Picard is a collection of scripts. This MultiQC module supports some but not all.
    The code for each script is split into its own file and adds a section to the
    module output if logs are found.
    """

    def __init__(
        self,
        name="Picard",
        anchor="picard",
        href="http://broadinstitute.github.io/picard/",
        info="is a set of Java command line tools for manipulating high-throughput sequencing data.",
        # No DOI to cite // doi=
        tools=tuple(TOOLS),
    ):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name=name,
            anchor=anchor,
            href=href,
            info=info,
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = dict()
        self.general_stats_data = dict()
        self.samples_parsed_by_tool = dict()

        for tool in tools:
            log.debug(f"Running picard tool {tool}")
            mod = globals()[tool]
            func = getattr(mod, "parse_reports", None)
            if func is not None:
                self.samples_parsed_by_tool[tool] = func(self)
                if len(self.samples_parsed_by_tool[tool]) > 0:
                    log.info(f"Found {len(self.samples_parsed_by_tool[tool])} {tool} reports")

        # Exit if we didn't find anything
        if all(len(v) == 0 for v in self.samples_parsed_by_tool.values()):
            raise ModuleNoSamplesFound
