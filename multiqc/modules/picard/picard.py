""" MultiQC module to parse output from Picard """


import logging
import os
import re
from collections import OrderedDict
from typing import Dict, Optional

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

        for tool in self.available_tools():
            func = getattr(globals()[tool], "parse_reports", None)
            if func is not None:
                n[tool] = func(self)
                if n[tool] > 0:
                    log.info(f"Found {n[tool]} {tool} reports")

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            raise ModuleNoSamplesFound

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)

    def available_tools(self):
        return [
            "AlignmentSummaryMetrics",
            "BaseDistributionByCycleMetrics",
            "CrosscheckFingerprints",
            "GcBiasMetrics",
            "HsMetrics",
            "InsertSizeMetrics",
            "MarkDuplicates",
            "OxoGMetrics",
            "QualityByCycleMetrics",
            "QualityScoreDistributionMetrics",
            "QualityYieldMetrics",
            "RnaSeqMetrics",
            "RrbsSummaryMetrics",
            "TargetedPcrMetrics",
            "VariantCallingMetrics",
            "ValidateSamFile",
            "WgsMetrics",
            "CollectIlluminaBasecallingMetrics",
            "CollectIlluminaLaneMetrics",
            "ExtractIlluminaBarcodes",
            "MarkIlluminaAdapters",
        ]

    def is_line_right_before_table(self, line: str) -> bool:
        """
        Picard logs from different samples can be concatenated together, so the module
        needs to know a marker to find where new sample information starts.

        Many tools and platforms - e.g. Sentieon, supported by MultiQC - use Picard
        internally for QC, while adding their own headers. So we can't assume that
        the headers will be the same for the same, and need to allow to override
        this function.
        """
        return line.startswith("## METRICS CLASS")

    def extract_sample_name(self, line: str, f: Dict) -> Optional[str]:
        """
        A file can be concatenated from multiple samples, so we can't just extract
        the sample name from the file name, and need a way to find sample name in the
        header. The copy of the originally used command is the best bet, as it's
        usually logged by Picard.
        """
        if line.startswith("# picard.analysis.") and ("INPUT=" in line or "INPUT" in line.split()):
            # Pull sample name from the input file name, recorded in the command line:
            fn_search = re.search(r"INPUT(?:=|\s+)(\[?[^\s]+\]?)", line, flags=re.IGNORECASE)
            if fn_search:
                s_name = os.path.basename(fn_search.group(1).strip("[]"))
                s_name = self.clean_s_name(s_name, f)
                return s_name
        return None

    @staticmethod
    def multiply_hundred(val):
        try:
            val = float(val) * 100
        except ValueError:
            pass
        return val
