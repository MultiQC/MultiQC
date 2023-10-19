""" MultiQC module to parse output from Picard """


import logging
import os
import re
from collections import OrderedDict
from typing import Dict, List, Optional, Union

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

        for tool, mod in self.available_tools().items():
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

    def available_tools(self) -> Dict:
        return {
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
        }

    @staticmethod
    def is_line_right_before_table(
        line: str,
        picard_class: Union[str, List[str]],
        sentieon_algo: Optional[str] = None,
    ) -> bool:
        """
        Picard logs from different samples can be concatenated together, so the module
        needs to know a marker to find where new sample information starts.

        The command line Picard tools themselves use the "## METRICS CLASS" header
        line for that purpose; however, the Picard classes often used by other
        tools and platforms - e.g. Sentieon and Parabricks  - while adding their own
        headers, so we need to handle them as well.
        """
        picard_classes = [picard_class] if isinstance(picard_class, str) else picard_class
        return (
            line.startswith("## METRICS CLASS")
            and any(c in line for c in picard_classes)
            or sentieon_algo
            and line.startswith("#SentieonCommandLine:")
            and f" --algo {sentieon_algo}" in line
        )

    def extract_sample_name(
        self,
        line: str,
        f: Dict,
        picard_tool: str,
        sentieon_algo: Optional[str] = None,
    ) -> Optional[str]:
        """
        Picard logs from different samples can be concatenated together, so we can't
        just extract the sample name from the file name. We need an alternative way
        to find sample names in the header lines instead. Picard records the command
        originally used to invoke itself, which is the best bet. Sentieon does the same,
        but slightly differently.
        """
        picard_command = line.startswith("# ") and ("INPUT=" in line or "INPUT" in line.split()) and picard_tool in line
        sentieon_command = (
            sentieon_algo
            and line.startswith("#SentieonCommandLine:")
            and f" --algo {sentieon_algo}" in line
            and " -i " in line
        )
        # Pull sample name from the input file name, recorded in the command line:
        fn_search = None
        if picard_command:
            fn_search = re.search(r"INPUT(?:=|\s+)(\[?[^\s]+\]?)", line, flags=re.IGNORECASE)
        elif sentieon_command:
            fn_search = re.search(r" -i\s+(\[?\S+\]?)", line, flags=re.IGNORECASE)
        if fn_search:
            f_name = os.path.basename(fn_search.group(1).strip("[]"))
            s_name = self.clean_s_name(f_name, f)
            return s_name
        return None

    @staticmethod
    def multiply_hundred(val):
        try:
            val = float(val) * 100
        except ValueError:
            pass
        return val
