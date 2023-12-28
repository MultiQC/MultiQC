# coding: utf-8
""" MultiQC submodule to parse output from Samtools markdup """

import logging
import re

from typing import Any
from typing import Callable
from typing import Dict
from typing import Tuple
from typing import Union

from multiqc.plots import bargraph
from multiqc.plots import table
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.modules.base_module import ModuleNoSamplesFound


# Initialise the logger
log = logging.getLogger(__name__)

class MarkdupReportMixin:
    """A MultiQC module for `samtools markdup` metrics."""

    search_key: str = "samtools/markdup"
    """The configuration key for storing search patterns about `samtools markdup` outputs."""

    @staticmethod
    def markdup_metric_patterns() -> Dict[str, Tuple[re.Pattern, Callable[[Any], Union[int, str]]]]:
        """Patterns for parsing the metrics within `samtools markdup` outputs."""
        return {
            "command": (re.compile(r"COMMAND:\s(.*)"), str),
            "optical_duplicate_distance": (re.compile(r"COMMAND:\s.*-d\s(\d+)"), int),
            "read": (re.compile(r"READ:\s(\d+)"), int),
            "written": (re.compile(r"WRITTEN:\s(\d+)"), int),
            "excluded": (re.compile(r"EXCLUDED:\s(\d+)"), int),
            "examined": (re.compile(r"EXAMINED:\s(\d+)"), int),
            "paired": (re.compile(r"PAIRED:\s(\d+)"), int),
            "single": (re.compile(r"SINGLE:\s(\d+)"), int),
            "duplicate_pair": (re.compile(r"DUPLICATE\sPAIR:\s(\d+)"), int),
            "duplicate_single": (re.compile(r"DUPLICATE\sSINGLE:\s(\d+)"), int),
            "duplicate_pair_optical": (re.compile(r"DUPLICATE\sPAIR\sOPTICAL:\s(\d+)"), int),
            "duplicate_single_optical": (re.compile(r"DUPLICATE\sSINGLE\sOPTICAL:\s(\d+)"), int),
            "duplicate_non_primary": (re.compile(r"DUPLICATE\sNON\sPRIMARY:\s(\d+)"), int),
            "duplicate_non_primary_optical": (re.compile(r"DUPLICATE\sNON\sPRIMARY\sOPTICAL:\s(\d+)"), int),
            "duplicate_primary_total": (re.compile(r"DUPLICATE\sPRIMARY\sTOTAL:\s(\d+)"), int),
            "duplicate_total": (re.compile(r"DUPLICATE\sTOTAL:\s(\d+)"), int),
            "estimated_library_size": (re.compile(r"ESTIMATED_LIBRARY_SIZE:\s(\d+)"), int),
        }

    @staticmethod
    def parse_contents(contents: str) -> Dict[str, Union[int, str]]:
        """Parse the contents of a `samtools markdup` output file."""
        metrics: Dict[str, Union[int, str]] = dict()

        for metric, (regex, converter) in MarkdupReportMixin.markdup_metric_patterns().items():
            match = regex.search(contents)
            if match:
                metrics[metric] = converter(match.group(1))

        return metrics

    @staticmethod
    def derive_metrics(metrics: Dict[str, Union[int, str]]) -> Dict[str, Union[int, str]]:
        """Derive custom metrics from the contents of a `samtools markdup` output file."""
        reads: int = metrics["paired"] + metrics["single"]
        metrics["duplicate_optical_total"] = metrics["duplicate_pair_optical"] + metrics["duplicate_single_optical"]
        metrics["duplicate_optical_fraction"] = reads and metrics["duplicate_optical_total"] / reads or 0.0
        metrics["duplicate_fraction"] = reads and metrics["duplicate_total"] / reads or 0.0
        metrics["duplicate_paired_non_optical"] = metrics["duplicate_pair"] - metrics["duplicate_pair_optical"]
        metrics["duplicate_single_non_optical"] = metrics["duplicate_single"] - metrics["duplicate_single_optical"]
        metrics["duplicate_non_primary_non_optical"] = metrics["duplicate_non_primary"] - metrics["duplicate_non_primary_optical"]
        metrics["non_duplicate"] = metrics["paired"] + metrics["single"] - metrics["duplicate_total"]
        return metrics

    def parse_samtools_markdup(self):
        """Initialize the MultiQC `samtools markdup` module."""
        super(MarkdupReportMixin, self).__init__(
            name="Samtools (Custom)",
            anchor="Samtools",
            target="Samtools",
            href="http://www.htslib.org",
            info=" is a suite of programs for interacting with high-throughput sequencing data.",
            doi="10.1093/bioinformatics/btp352",
        )

        metrics: Dict[str, Dict[str, Union[int, str]]] = dict()

        # TODO: Check support for the JSON output form of this log file
        # TODO: Support parsing from a stdout log file in either text/JSON mode
        for file in self.find_log_files(self.search_key):
            parsed = self.parse_contents(file["f"])
            metrics[file["s_name"]] = self.derive_metrics(parsed)

        metrics = self.ignore_samples(data=metrics)

        if len(metrics) == 0:
            return 0

        self.add_software_version(None)

        self.write_data_file(data=metrics, fn="multiqc_samtools")

        # General Statistics ###########################################################

        headers = {
            "duplicate_fraction": {
                "title": "% Duplicates",
                "description": "The percent of all types of duplicate reads.",
                "min": 0,
                "modify": lambda x: x * 100,
                "format": "{:,.0f}",
                "suffix": "%",
                "scale": "RdYlGn-rev",
            },
            "estimated_library_size": {
                "title": "Estimated Library Size",
                "description": "The estimated library size after de-duplication.",
                "min": 0,
                "format": "{:,.0f}",
            },
        }

        self.general_stats_addcols(data=metrics, headers=headers)

        # Detailed Metrics #############################################################

        headers = {
            "estimated_library_size": {
                "title": "Estimated Library Size",
                "description": "The estimated library size after de-duplication.",
                "min": 0,
                "format": "{:,.0f}",
            },
            "optical_duplicate_distance": {
                "title": "Optical Distance",
                "description": "The optical distance for considering instrument duplicates.",
                "min": 0,
                "format": "{:,.0f}",
            },
            "duplicate_fraction": {
                "title": "% Duplicates",
                "description": "The percent of all types of duplicate reads.",
                "min": 0,
                "modify": lambda x: x * 100,
                "format": "{:,.0f}",
                "suffix": "%",
                "scale": "RdYlGn-rev",
            },
            "duplicate_optical_fraction": {
                "title": "% Optical Duplicates",
                "description": "The percent of optical/clustering duplicate reads.",
                "min": 0,
                "modify": lambda x: x * 100,
                "format": "{:,.0f}",
                "suffix": "%",
                "scale": "RdYlGn-rev",
            },
        }

        self.add_section(
            name="Samtools markdup",
            anchor="samtools-markdup",
            description=(
                "Optical duplicates are due to either optical or clustering-based artifacts. "
                + "See the following links to learn more about instrument-based duplicate "
                + "artifacts:"
                + "<br>"
                + "<ul>"
                + '<li><a href="https://core-genomics.blogspot.com/2016/05/increased-read-duplication-on-patterned.html" '
                + 'target="_blank">Core Genomics Post: Increased Read Duplication on Patterned '
                + "Flowcells</a>"
                + "</li>"
                + '<li><a href="https://sequencing.qcfail.com/articles/illumina-patterned-flow-cells-generate-duplicated-sequences/" '
                + 'target="_blank">QC Fail Post: Illumina Patterned Flow Cells Generate '
                + "Duplicated Sequences</a>"
                + "</li>"
                + "</ul>."
            ),
            plot=table.plot(
                data=metrics,
                headers=headers,
                pconfig={
                    "namespace": self.name,
                    "id": "samtools-markdup-table",
                    "table_title": "Samtools: Duplicate Marked SAM Records (Alignments)",
                    "sortRows": False,
                },
            ),
        )

        # Stacked Bar Plot #############################################################

        bargraph_config = {
            "namespace": self.name,
            "id": "samtools-markdup-fraction",
            "title": "Samtools: Duplicate Marked SAM Records (Alignments)",
            "cpswitch": True,
            "ylab": "% SAM Records",
        }

        # TODO; Confirm that under all ways to run samtools markdup, these sum up to "all records"
        keys: Dict[str, Dict[str, str]] = {
            "non_duplicate": {"name": "Non-Duplicates"},
            "duplicate_pair_optical": {"name": "Optical Duplicates in Pairs"},
            "duplicate_single_optical": {"name": "Optical Duplicates in Singletons"},
            "duplicate_non_primary_optical": {"name": "Optical Non-Primary Duplicate"},
            "duplicate_pair_non_optical": {"name": "Non-optical Duplicates in Pairs"},
            "duplicate_single_non_optical": {"name": "Non-optical Duplicates in Singletons"},
            "duplicate_non_primary_non_optical": {"name": "Non-Optical Non-Primary Duplicates"},
            "excluded": {"name": "Ignored (QC fail or unmapped)"},
        }

        self.add_section(
            description=(
                "For more information about the duplicate categories, see the "
                + '<a href="https://www.htslib.org/doc/samtools-markdup.html#STATISTICS" '
                + 'target="_blank">samtools documentation</a>. '
            ),
            plot=bargraph.plot(data=metrics, cats=keys, pconfig=bargraph_config),
        )

        return len(metrics)
