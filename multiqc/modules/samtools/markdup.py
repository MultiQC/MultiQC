# coding: utf-8
""" MultiQC submodule to parse output from Samtools markdup """

import logging
import re

from typing import Callable
from typing import Dict
from typing import Tuple
from typing import Union

from multiqc.plots import bargraph
from multiqc.plots import table

log = logging.getLogger(__name__)


class MarkdupReportMixin:
    """A MultiQC module for `samtools markdup` metrics."""

    @staticmethod
    def markdup_metric_patterns() -> Dict[str, Tuple[re.Pattern, Callable[[str], Union[int, str]]]]:
        """Patterns for parsing the metrics within `samtools markdup` outputs."""
        return {
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
    def parse_contents(contents: str) -> Dict[str, Union[int, float]]:
        """
        Parse the contents of a `samtools markdup` output file.
        """
        data: Dict[str, Union[int, float]] = dict()
        for metric, (regex, converter) in MarkdupReportMixin.markdup_metric_patterns().items():
            match = regex.search(contents)
            if match:
                data[metric] = converter(match.group(1))

        return data

    def parse_samtools_markdup(self):
        val_by_metric_by_sample: Dict[str, Dict[str, Union[int, float]]] = dict()

        for f in self.find_log_files("samtools/markdup"):
            data = self.parse_contents(f["f"])
            if not data:
                continue
            if f["s_name"] in val_by_metric_by_sample:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {f['s_name']}")
            self.add_data_source(f, section="markdup")

            # Derive more metrics from counts
            n_reads: int = data["paired"] + data["single"]
            data["duplicate_optical_total"] = data["duplicate_pair_optical"] + data["duplicate_single_optical"]
            data["duplicate_optical_fraction"] = data["duplicate_optical_total"] / n_reads if n_reads > 0 else 0.0
            data["duplicate_fraction"] = data["duplicate_total"] / n_reads if n_reads > 0 else 0.0
            data["duplicate_paired_non_optical"] = data["duplicate_pair"] - data["duplicate_pair_optical"]
            data["duplicate_single_non_optical"] = data["duplicate_single"] - data["duplicate_single_optical"]
            data["duplicate_non_primary_non_optical"] = (
                data["duplicate_non_primary"] - data["duplicate_non_primary_optical"]
            )
            data["non_duplicate"] = data["paired"] + data["single"] - data["duplicate_total"]
            val_by_metric_by_sample[f["s_name"]] = data

        val_by_metric_by_sample = self.ignore_samples(val_by_metric_by_sample)
        if len(val_by_metric_by_sample) == 0:
            return 0

        # Superfluous function call to confirm that it is used in this module
        self.add_software_version(None)

        self.write_data_file(val_by_metric_by_sample, fn="multiqc_samtools_markdup")

        genstats_headers = {
            "duplicate_fraction": {
                "title": "Duplicates",
                "description": "The percent of all types of duplicate reads",
                "min": 0,
                "max": 100,
                "modify": lambda x: x * 100,
                "suffix": "%",
                "scale": "OrRd",
            },
            "estimated_library_size": {
                "title": "Est. library size",
                "description": "The estimated library size after de-duplication.",
                "min": 0,
                "format": "{:,d}",
            },
        }
        self.general_stats_addcols(data=val_by_metric_by_sample, headers=genstats_headers, namespace="markdup")

        self.add_section(
            name="Samtools markdup: stats",
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
                data=val_by_metric_by_sample,
                headers=dict(
                    genstats_headers,
                    **{
                        "optical_duplicate_distance": {
                            "title": "Optical distance",
                            "description": "The optical distance for considering instrument duplicates",
                            "min": 0,
                            "format": "{:,d}",
                            "scale": "RdYlGn",
                        },
                        "duplicate_optical_fraction": {
                            "title": "Optical dups",
                            "description": "The percent of optical/clustering duplicate reads",
                            "min": 0,
                            "max": 100,
                            "modify": lambda x: x * 100,
                            "suffix": "%",
                            "scale": "RdYlGn-rev",
                        },
                    },
                ),
                pconfig={
                    "id": "samtools-markdup-table",
                    "title": "Samtools: duplicate-marked SAM records (alignments)",
                },
            ),
        )

        # Stacked Bar Plot
        pconfig = {
            "id": "samtools-markdup-fraction",
            "title": "Samtools markdup: duplicate categories",
            "ylab": "SAM Records",
        }

        keys: Dict[str, Dict[str, str]] = {
            "non_duplicate": {"name": "Non-duplicates"},
            "duplicate_pair_optical": {"name": "Optical duplicates in pairs"},
            "duplicate_single_optical": {"name": "Optical duplicates in singletons"},
            "duplicate_non_primary_optical": {"name": "Optical non-primary duplicate"},
            "duplicate_pair_non_optical": {"name": "Non-optical duplicates in pairs"},
            "duplicate_single_non_optical": {"name": "Non-optical duplicates in singletons"},
            "duplicate_non_primary_non_optical": {"name": "Non-optical non-primary duplicates"},
            "excluded": {"name": "Ignored (QC fail or unmapped)"},
        }

        self.add_section(
            name="Samtools markdup: duplicate categories",
            anchor="samtools-markdup-categories",
            description=(
                "For more information about the duplicate categories, see the "
                + '<a href="https://www.htslib.org/doc/samtools-markdup.html#STATISTICS" '
                + 'target="_blank">samtools documentation</a>. '
            ),
            plot=bargraph.plot(data=val_by_metric_by_sample, cats=keys, pconfig=pconfig),
        )

        return len(val_by_metric_by_sample)
