import csv
import logging
from typing import Any, TextIO, cast

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table
from multiqc.plots.table_object import ColumnDict

log = logging.getLogger(__name__)

INTEGER_FIELDS = {
    "ambiguous_reads",
    "assigned_reads",
    "exact_reads",
    "invalid_reads",
    "k1_rescued_reads",
    "no_match_reads",
    "targets_observed",
    "zero_count_targets",
    "candidates_verified",
    "total_count",
}

FLOAT_FIELDS = {
    "assignment_rate",
    "exact_rate",
    "rescue_rate",
    "ambiguous_rate",
    "no_match_rate",
    "coverage_fraction",
    "zero_count_fraction",
    "gini_index",
    "top_1pct_fraction",
}


def _coerce(field: str, value: str) -> Any:
    """Convert the numeric fields used by the DotMatch report schemas."""
    if value == "":
        return ""
    if field in INTEGER_FIELDS:
        try:
            return int(value)
        except ValueError:
            return value
    if field in FLOAT_FIELDS:
        try:
            return float(value)
        except ValueError:
            return value
    return value


def _read_tsv(handle: TextIO) -> list[dict[str, Any]]:
    """Read a DotMatch TSV report and coerce its known numeric columns."""
    return [
        {key: _coerce(key, value or "") for key, value in row.items()}
        for row in csv.DictReader(handle, delimiter="\t")
        if any(value not in (None, "") for value in row.values())
    ]


def _parse_sample_qc(handle: TextIO, fallback_sample: str) -> dict[str, dict[str, Any]]:
    data: dict[str, dict[str, Any]] = {}
    for row in _read_tsv(handle):
        sample = str(row.get("sample_id") or fallback_sample)
        data[sample] = {
            "total_reads": row.get("total_reads", ""),
            "valid_extracted_reads": row.get("valid_extracted_reads", ""),
            "assigned_reads": row.get("assigned_reads", ""),
            "exact_reads": row.get("exact_reads", ""),
            "k1_rescued_reads": row.get("k1_rescued_reads", ""),
            "ambiguous_reads": row.get("ambiguous_reads", ""),
            "no_match_reads": row.get("no_match_reads", ""),
            "invalid_reads": row.get("invalid_reads", ""),
            "assignment_rate": row.get("assignment_rate", ""),
            "exact_rate": row.get("exact_rate", ""),
            "rescue_rate": row.get("rescue_rate", ""),
            "ambiguous_rate": row.get("ambiguous_rate", ""),
            "no_match_rate": row.get("no_match_rate", ""),
            "targets_observed": row.get("targets_observed", ""),
            "zero_count_targets": row.get("zero_count_targets", ""),
            "gini_index": row.get("gini_index", ""),
            "top_1pct_read_fraction": row.get("top_1pct_read_fraction", ""),
            "candidates_verified": row.get("candidates_verified", ""),
        }
    return data


def _parse_crispr_qc(handle: TextIO, fallback_sample: str) -> dict[str, dict[str, Any]]:
    data: dict[str, dict[str, Any]] = {}
    for row in _read_tsv(handle):
        sample = str(row.get("sample_id") or fallback_sample)
        data[sample] = {
            "qc_status": row.get("qc_status", ""),
            "total_count": row.get("total_count", ""),
            "coverage_fraction": row.get("coverage_fraction", ""),
            "zero_count_fraction": row.get("zero_count_fraction", ""),
            "gini_index": row.get("gini_index", ""),
            "top_1pct_fraction": row.get("top_1pct_fraction", ""),
            "assignment_rate": row.get("assignment_rate", ""),
            "ambiguous_rate": row.get("ambiguous_rate", ""),
            "no_match_rate": row.get("no_match_rate", ""),
            "invalid_rate": row.get("invalid_rate", ""),
        }
    return data


class MultiqcModule(BaseMultiqcModule):
    """
    Summarise DotMatch assignment and CRISPR guide-library QC reports.

    The module recognises the stable ``sample_qc.tsv`` and
    ``crispr_qc.summary.tsv`` artifacts written by DotMatch. Ambiguous reads
    remain a separate metric and are never folded into assigned counts.
    """

    def __init__(self):
        super().__init__(
            name="DotMatch",
            anchor="dotmatch",
            href="https://github.com/dnncha/dotmatch",
            info=("Deterministic short-DNA known-target assignment with explicit ambiguity handling."),
        )

        self.sample_qc: dict[str, dict[str, Any]] = {}
        for file_data in self.find_log_files("dotmatch/sample_qc", filecontents=False, filehandles=True):
            self.sample_qc.update(_parse_sample_qc(cast(TextIO, file_data["f"]), file_data["s_name"]))
            self.add_data_source(file_data, file_data["s_name"])

        self.crispr_qc: dict[str, dict[str, Any]] = {}
        for file_data in self.find_log_files("dotmatch/crispr_qc", filecontents=False, filehandles=True):
            self.crispr_qc.update(_parse_crispr_qc(cast(TextIO, file_data["f"]), file_data["s_name"]))
            self.add_data_source(file_data, file_data["s_name"])

        self.sample_qc = self.ignore_samples(self.sample_qc)
        self.crispr_qc = self.ignore_samples(self.crispr_qc)

        if not self.sample_qc and not self.crispr_qc:
            raise ModuleNoSamplesFound

        if self.sample_qc:
            log.info("Found %d DotMatch sample QC rows", len(self.sample_qc))
            sample_headers: dict[str, ColumnDict] = {
                "assigned_reads": {
                    "title": "Assigned",
                    "description": "Reads assigned to exactly one known target.",
                    "format": "{:,.0f}",
                    "scale": "Blues",
                },
                "exact_reads": {"title": "Exact", "format": "{:,.0f}", "scale": "Blues"},
                "k1_rescued_reads": {
                    "title": "Rescued",
                    "format": "{:,.0f}",
                    "scale": "Oranges",
                },
                "ambiguous_reads": {
                    "title": "Ambiguous",
                    "description": "Reads within the configured radius of multiple targets.",
                    "format": "{:,.0f}",
                    "scale": "Reds",
                },
                "no_match_reads": {"title": "No match", "format": "{:,.0f}", "scale": "Greys"},
                "assignment_rate": {
                    "title": "Assign %",
                    "format": "{:.2%}",
                    "min": 0,
                    "max": 1,
                    "scale": "YlGn",
                },
                "ambiguous_rate": {
                    "title": "Ambig %",
                    "format": "{:.2%}",
                    "min": 0,
                    "max": 1,
                    "scale": "OrRd",
                },
                "no_match_rate": {
                    "title": "No match %",
                    "format": "{:.2%}",
                    "min": 0,
                    "max": 1,
                    "scale": "Greys",
                },
            }
            self.general_stats_addcols(
                self.sample_qc,
                {
                    "assignment_rate": {
                        "title": "Assign %",
                        "description": "Fraction of valid reads assigned to one known target.",
                        "min": 0,
                        "max": 1,
                        "format": "{:.2%}",
                        "scale": "YlGn",
                    },
                    "ambiguous_rate": {
                        "title": "Ambig %",
                        "description": "Fraction of reads excluded because multiple targets were within the radius.",
                        "min": 0,
                        "max": 1,
                        "format": "{:.2%}",
                        "scale": "OrRd",
                    },
                },
            )
            self.add_section(
                name="DotMatch Sample QC",
                anchor="dotmatch_sample_qc",
                description=(
                    "Per-sample assignment outcomes. Ambiguous reads remain separate from uniquely assigned reads."
                ),
                plot=table.plot(
                    self.sample_qc,
                    sample_headers,
                    pconfig={"id": "dotmatch_sample_qc_table", "title": "DotMatch Sample QC"},
                ),
            )
            self.write_data_file(self.sample_qc, "multiqc_dotmatch_sample_qc")

        if self.crispr_qc:
            log.info("Found %d DotMatch CRISPR QC rows", len(self.crispr_qc))
            self.add_section(
                name="DotMatch CRISPR QC",
                anchor="dotmatch_crispr_qc",
                description=("Guide-library representation and assignment QC from DotMatch CRISPR count reports."),
                plot=table.plot(
                    self.crispr_qc,
                    cast(
                        dict[str, ColumnDict],
                        {
                            "qc_status": {"title": "QC"},
                            "total_count": {"title": "Total guides", "format": "{:,.0f}", "scale": "Blues"},
                            "coverage_fraction": {
                                "title": "Coverage",
                                "format": "{:.2%}",
                                "min": 0,
                                "max": 1,
                                "scale": "YlGn",
                            },
                            "zero_count_fraction": {
                                "title": "Zero %",
                                "format": "{:.2%}",
                                "min": 0,
                                "max": 1,
                                "scale": "OrRd",
                            },
                            "gini_index": {"title": "Gini", "format": "{:.3f}", "scale": "Purples"},
                            "assignment_rate": {
                                "title": "Assign %",
                                "format": "{:.2%}",
                                "min": 0,
                                "max": 1,
                                "scale": "YlGn",
                            },
                        },
                    ),
                    pconfig={"id": "dotmatch_crispr_qc_table", "title": "DotMatch CRISPR QC"},
                ),
            )
            self.write_data_file(self.crispr_qc, "multiqc_dotmatch_crispr_qc")
