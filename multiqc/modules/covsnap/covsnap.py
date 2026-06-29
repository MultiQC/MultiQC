import json
import logging
from typing import Dict, Optional, Union

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)

# Classification status severity, from least to most severe. Used to pick the
# worst status across a sample's targets and to order the stacked bar plot.
STATUS_SEVERITY = {
    "PASS": 0,
    "LOW_COVERAGE": 1,
    "LOW_EXON": 2,
    "UNEVEN": 3,
    "DROP_OUT": 4,
}

STATUS_COLORS = {
    "PASS": "#5cb85c",
    "LOW_COVERAGE": "#f0ad4e",
    "LOW_EXON": "#fee08b",
    "UNEVEN": "#5bc0de",
    "DROP_OUT": "#d9534f",
}


class MultiqcModule(BaseMultiqcModule):
    """
    The module parses the JSON report produced by covsnap with `--format json`
    (a file named `<output>.json`). Each file is recognised by its
    `covsnap_version` key, and the sample name is read from the `run.sample_name`
    field inside the JSON rather than from the file name.

    covsnap classifies each target (gene, region, or BED interval) as one of
    PASS, LOW_COVERAGE, UNEVEN, LOW_EXON, or DROP_OUT. The module reports, per
    sample:

    - the fraction of targets that PASS,
    - the length-weighted mean depth,
    - the length-weighted percentage of bases at >=20x,
    - the worst (most severe) status across all targets.

    Two sections are added: a per-target coverage table (one row per
    sample/target) and a stacked bar plot of target counts by classification
    status. To generate the input, run covsnap with the JSON output enabled, for
    example:

    ```
    covsnap sample.bam BRCA1,TP53 --format json
    ```
    """

    def __init__(self):
        super().__init__(
            name="covsnap",
            anchor="covsnap",
            href="https://github.com/enes-ak/covsnap",
            info="Computes per-target and per-exon depth-of-coverage metrics for targeted sequencing QC, with automated PASS/FAIL classification",
            doi="10.5281/zenodo.18732742",
        )

        # Per-sample summary metrics (general stats + summary table)
        summary_by_sample: Dict[str, Dict[str, Union[float, int, str]]] = dict()
        # Per-target metrics, keyed "sample: target" (detail table)
        target_rows: Dict[str, Dict[str, Union[float, int, str]]] = dict()
        # Per-sample counts of targets by classification status (bar plot)
        status_counts_by_sample: Dict[str, Dict[str, int]] = dict()

        for f in self.find_log_files("covsnap", filehandles=False):
            parsed = self.parse_covsnap_json(f)
            if parsed is None:
                continue
            s_name, summary, rows, status_counts = parsed
            if s_name in summary_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            summary_by_sample[s_name] = summary
            status_counts_by_sample[s_name] = status_counts
            for row_key, row in rows.items():
                target_rows[row_key] = row
            self.add_data_source(f, s_name)

        summary_by_sample = self.ignore_samples(summary_by_sample)
        if len(summary_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(summary_by_sample)} reports")

        self.general_stats_table(summary_by_sample)
        self.target_table(target_rows)
        self.status_bar_plot(status_counts_by_sample)

        self.write_data_file(summary_by_sample, "multiqc_covsnap")

    def parse_covsnap_json(self, f):
        """
        Parse one covsnap JSON report. Returns a tuple of
        (sample_name, summary_dict, per_target_rows, status_counts) or None if
        the file is not a covsnap JSON report.
        """
        try:
            content = json.loads(f["f"])
        except (json.JSONDecodeError, TypeError):
            return None
        if not isinstance(content, dict) or "covsnap_version" not in content:
            return None

        run = content["run"]
        summary = content["summary"]
        targets = content["targets"]

        s_name = self.clean_s_name(run["sample_name"], f)
        self.add_software_version(content["covsnap_version"], sample=s_name)

        n_targets = summary["n_targets"]
        n_pass = summary["n_pass"]

        # Length-weighted aggregates across targets.
        total_len = sum(t["length_bp"] for t in targets) or 1
        mean_depth = sum(t["mean_depth"] * t["length_bp"] for t in targets) / total_len

        summary_row: Dict[str, Union[float, int, str]] = {
            "pct_targets_pass": (n_pass / n_targets * 100) if n_targets else 0.0,
            "mean_depth": round(mean_depth, 2),
            "n_targets": n_targets,
            "n_pass": n_pass,
            "worst_status": _worst_status(summary["status_counts"]),
        }
        # pct_thresholds[20] is optional (only present when 20 is a configured threshold).
        if any("20" in t["pct_thresholds"] for t in targets):
            pct20 = sum(t["pct_thresholds"].get("20", 0.0) * t["length_bp"] for t in targets) / total_len
            summary_row["pct_ge_20"] = round(pct20, 2)

        rows: Dict[str, Dict[str, Union[float, int, str]]] = dict()
        for t in targets:
            row: Dict[str, Union[float, int, str]] = {
                "contig": t["contig"],
                "length_bp": t["length_bp"],
                "mean_depth": t["mean_depth"],
                "median_depth": t["median_depth"],
                "pct_zero": t["pct_zero"],
                "coverage_status": t["coverage_status"],
            }
            if "20" in t["pct_thresholds"]:
                row["pct_ge_20"] = t["pct_thresholds"]["20"]
            rows[f"{s_name}: {t['target_id']}"] = row

        status_counts: Dict[str, int] = {status: 0 for status in STATUS_SEVERITY}
        for status, count in summary["status_counts"].items():
            status_counts[status] = status_counts.get(status, 0) + count

        return s_name, summary_row, rows, status_counts

    def general_stats_table(self, summary_by_sample):
        """Add covsnap summary metrics to the General Statistics table."""
        headers: Dict[str, Dict] = {
            "pct_targets_pass": {
                "title": "Targets PASS",
                "description": "Percentage of targets classified PASS by covsnap",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.1f}",
            },
            "pct_ge_20": {
                "title": "% >= 20x",
                "description": "Length-weighted percentage of target bases covered at >=20x",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "YlGn",
                "format": "{:,.1f}",
            },
            "mean_depth": {
                "title": "Mean depth",
                "description": "Length-weighted mean depth across targets",
                "min": 0,
                "suffix": "x",
                "scale": "BuPu",
                "format": "{:,.1f}",
            },
            "worst_status": {
                "title": "Worst status",
                "description": "Most severe covsnap classification across targets",
                "scale": False,
                "hidden": True,
            },
        }
        self.general_stats_addcols(summary_by_sample, headers)

    def target_table(self, target_rows):
        """Add a per-target coverage table (one row per sample/target)."""
        headers: Dict[str, Dict] = {
            "contig": {
                "title": "Contig",
                "description": "Target contig",
                "scale": False,
            },
            "mean_depth": {
                "title": "Mean depth",
                "description": "Mean depth across the target",
                "min": 0,
                "suffix": "x",
                "scale": "BuPu",
                "format": "{:,.1f}",
            },
            "median_depth": {
                "title": "Median depth",
                "description": "Median depth across the target",
                "min": 0,
                "suffix": "x",
                "scale": "BuPu",
                "format": "{:,.1f}",
            },
            "pct_ge_20": {
                "title": "% >= 20x",
                "description": "Percentage of target bases covered at >=20x",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "YlGn",
                "format": "{:,.1f}",
            },
            "pct_zero": {
                "title": "% zero",
                "description": "Percentage of target bases with zero coverage",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
                "format": "{:,.2f}",
            },
            "length_bp": {
                "title": "Length",
                "description": "Target length in base pairs",
                "min": 0,
                "suffix": " bp",
                "scale": "Greys",
                "format": "{:,d}",
                "hidden": True,
            },
            "coverage_status": {
                "title": "Status",
                "description": "covsnap classification for this target",
                "scale": False,
            },
        }
        self.add_section(
            name="Per-target coverage",
            anchor="covsnap-target-table",
            description="Depth-of-coverage metrics and PASS/FAIL classification for each target.",
            helptext="""
            Each row is one target (gene, region, or BED interval) within a sample.
            `% >= 20x` is the breadth of coverage at the 20x clinical threshold;
            `% zero` is the fraction of the target with no coverage. The `Status`
            column is covsnap's classification: PASS, LOW_COVERAGE, UNEVEN,
            LOW_EXON, or DROP_OUT.
            """,
            plot=table.plot(
                target_rows,
                headers,
                pconfig={
                    "id": "covsnap-target-table",
                    "title": "covsnap: Per-target coverage",
                    "col1_header": "Sample: target",
                },
            ),
        )

    def status_bar_plot(self, status_counts_by_sample):
        """Add a stacked bar plot of target counts by classification status."""
        cats = {
            status: {"name": status, "color": STATUS_COLORS[status]}
            for status in sorted(STATUS_SEVERITY, key=lambda s: STATUS_SEVERITY[s])
        }
        self.add_section(
            name="Target status distribution",
            anchor="covsnap-status-dist",
            description="Number of targets in each covsnap classification, per sample.",
            helptext="""
            covsnap classifies every target into one of five statuses. This plot
            shows how many targets fall into each status for every sample, making
            it easy to spot samples with many failing targets.

            * PASS: coverage is adequate across the target.
            * LOW_COVERAGE: breadth at 20x is below the configured threshold.
            * UNEVEN: high mean depth but highly variable distribution.
            * LOW_EXON: one or more exons have critically low coverage.
            * DROP_OUT: large portions of the target have no coverage.
            """,
            plot=bargraph.plot(
                status_counts_by_sample,
                cats,
                pconfig={
                    "id": "covsnap-status-dist",
                    "title": "covsnap: Target status distribution",
                    "ylab": "Number of targets",
                    "cpswitch_counts_label": "Number of targets",
                },
            ),
        )


def _worst_status(status_counts: Dict[str, int]) -> str:
    """Return the most severe status present in a status_counts mapping."""
    present = [s for s, n in status_counts.items() if n > 0]
    if not present:
        return ""
    return max(present, key=lambda s: STATUS_SEVERITY.get(s, 0))
