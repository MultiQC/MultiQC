""" MultiQC module to parse output from Picard """


import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Import the Picard submodules
from . import (
    AlignmentSummaryMetrics,
    BaseDistributionByCycleMetrics,
    CollectIlluminaBasecallingMetrics,
    CollectIlluminaLaneMetrics,
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

# Initialise the logger
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
    ):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name=name,
            anchor=anchor,
            href=href,
            info=info,
        )

        # Set up class objects to hold parsed data
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

        for tool, mod in {
            "AlignmentSummaryMetrics": AlignmentSummaryMetrics,
            "BaseDistributionByCycleMetrics": BaseDistributionByCycleMetrics,
            "CrosscheckFingerprints": CrosscheckFingerprints,
            "GcBiasMetrics": GcBiasMetrics,
            "HsMetrics": HsMetrics,
            "InsertSizeMetrics": InsertSizeMetrics,
            "MarkDuplicates": MarkDuplicates,
            "OxoGMetrics": OxoGMetrics,
            "QualityByCycleMetrics": QualityByCycleMetrics,
            "QualityScoreDistributionMetrics": QualityScoreDistributionMetrics,
            "QualityYieldMetrics": QualityYieldMetrics,
            "RnaSeqMetrics": RnaSeqMetrics,
            "RrbsSummaryMetrics": RrbsSummaryMetrics,
            "TargetedPcrMetrics": TargetedPcrMetrics,
            "VariantCallingMetrics": VariantCallingMetrics,
            "ValidateSamFile": ValidateSamFile,
            "WgsMetrics": WgsMetrics,
            "CollectIlluminaBasecallingMetrics": CollectIlluminaBasecallingMetrics,
            "CollectIlluminaLaneMetrics": CollectIlluminaLaneMetrics,
            "ExtractIlluminaBarcodes": ExtractIlluminaBarcodes,
            "MarkIlluminaAdapters": MarkIlluminaAdapters,
        }.items():
            func = getattr(mod, "parse_reports", None)
            if func is not None:
                n[tool] = func(self)
                if n[tool] > 0:
                    log.info(f"Found {n[tool]} {tool} reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)
