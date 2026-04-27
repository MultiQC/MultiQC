"""
Native MultiQC module for Trim Galore v2.x ("Oxidized Edition").

Parses the structured `*_trimming_report.json` file emitted by Trim Galore
v2.x alongside the legacy text report (schema v1).  See:
https://github.com/MultiQC/MultiQC/issues/3529

Why a dedicated module rather than reusing the `cutadapt` module:
  * Software-versions table shows "Trim Galore" with the correct version
    instead of the misleading "Cutadapt 4.0" backwards-compatibility shim.
  * Native access to TrimGalore-specific stats not present in Cutadapt
    output (RRBS truncation counts, poly-A / poly-G trimming summaries,
    paired-end pair-validation outcomes).
  * The JSON is purely additive — Trim Galore continues to emit the
    text report unchanged, so the existing `cutadapt` module path
    keeps working for users who don't enable this module.
"""

import json
import logging
import re
from typing import Any, Dict, Optional, Tuple

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph

log = logging.getLogger(__name__)


# ── Schema constants ────────────────────────────────────────────────────────
#
# Trim Galore JSON schema v1 — these keys are stable for the v2.x line.
# When the schema bumps, add handling for newer versions and keep these
# constants for the v1 path.

SCHEMA_VERSION_SUPPORTED = 1
TOOL_NAME = "Trim Galore"


# ── Module ──────────────────────────────────────────────────────────────────


class MultiqcModule(BaseMultiqcModule):
    """
    Trim Galore v2.x ("Oxidized Edition") — Rust rewrite of the original
    Perl script. Single binary, zero external runtime deps. Same CLI and
    same output filenames as v0.6.x; emits an additional structured JSON
    report (`*_trimming_report.json`) that this module parses.

    The legacy text report (`*_trimming_report.txt`) is still produced
    unchanged and continues to be parsed by the `cutadapt` module via the
    backwards-compatibility shim.  Users who want native Trim Galore
    parsing should disable the `cutadapt` module to avoid duplicate
    sample entries:

    ```yaml
    disable_modules:
      - cutadapt
    ```
    """

    def __init__(self):
        super().__init__(
            name="Trim Galore",
            anchor="trim_galore",
            href="https://github.com/FelixKrueger/TrimGalore",
            info=(
                "Quality and adapter trimming for next-generation sequencing data, "
                "with special handling for RRBS libraries."
            ),
            extra=(
                "Trim Galore v2.x is a Rust rewrite of the original Perl tool. "
                "It is a drop-in replacement for v0.6.x scripts and pipelines: "
                "same CLI, same output filenames, same text report. This MultiQC "
                "module parses the structured JSON report (`*_trimming_report.json`) "
                "introduced in v2.0 alongside the legacy text report."
            ),
            doi="10.5281/zenodo.5127899",
        )

        data_by_sample: Dict[str, Dict[str, Any]] = {}
        for f in self.find_log_files("trim_galore", filehandles=True):
            parsed = self._parse_log(f)
            if parsed is None:
                continue
            s_name, payload = parsed
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            data_by_sample[s_name] = payload
            self.add_data_source(f, s_name=s_name)

            # Software version — single canonical entry per sample.
            version = payload.get("trim_galore_version")
            if version:
                self.add_software_version(version, s_name)

        data_by_sample = self.ignore_samples(data_by_sample)
        if not data_by_sample:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(data_by_sample)} reports")

        # Persist the raw parsed JSON for downstream use / inspection.
        self.write_data_file(data_by_sample, "multiqc_trim_galore")

        # ── General stats table ──
        self._general_stats_table(data_by_sample)

        # ── Filtering disposition bar plot ──
        self.add_section(
            name="Filtered Reads",
            anchor="trim_galore_filtered_reads",
            description=(
                "Disposition of reads after trimming. Reads that passed all filters "
                "land in the trimmed output; the rest are categorised by which filter "
                "rejected them."
            ),
            plot=self._filtered_reads_plot(data_by_sample),
        )

        # ── Adapter length distribution ──
        # One trace per (sample, adapter) — collapses cleanly when most samples
        # use a single adapter.
        adapter_plot = self._adapter_length_plot(data_by_sample)
        if adapter_plot is not None:
            self.add_section(
                name="Adapter Length Distribution",
                anchor="trim_galore_adapter_lengths",
                description=(
                    "Histogram of adapter match lengths per read. Tall left tail at "
                    "1 bp is normal — Trim Galore's default `--stringency 1` accepts "
                    "single-base overlaps, which represent random hits rather than "
                    "real adapter contamination."
                ),
                plot=adapter_plot,
            )

    # ── Parsing ─────────────────────────────────────────────────────────

    def _parse_log(self, f) -> Optional[Tuple[str, Dict[str, Any]]]:
        """Parse a single TrimGalore JSON report. Returns (sample_name, payload)
        or None if the file is malformed / wrong schema."""
        try:
            payload = json.load(f["f"])
        except json.JSONDecodeError as e:
            log.warning(f"Could not parse {f['fn']!r}: {e}")
            return None

        if payload.get("tool") != TOOL_NAME:
            return None
        schema = payload.get("schema_version")
        if schema != SCHEMA_VERSION_SUPPORTED:
            log.warning(
                f"Skipping {f['fn']!r}: unsupported schema_version {schema!r} "
                f"(expected {SCHEMA_VERSION_SUPPORTED})"
            )
            return None

        # Sample name: PE reports list BOTH input filenames in
        # `input_filenames`, so deriving from index 0 collapses R1 + R2 onto
        # the same sample. Use `read_number` (1 or 2) to pick the matching
        # entry. SE reports have `read_number: 1` and `input_filenames` of
        # length 1, so the same logic gives the correct SE behaviour.
        input_filenames = payload.get("input_filenames") or []
        read_number = payload.get("read_number") or 1
        if not input_filenames:
            log.warning(f"Skipping {f['fn']!r}: no input_filenames")
            return None
        idx = max(0, min(read_number - 1, len(input_filenames) - 1))
        s_name = self.clean_s_name(_strip_fastq_suffix(input_filenames[idx]), f)
        return s_name, payload

    # ── Section: general stats ──────────────────────────────────────────

    def _general_stats_table(self, data_by_sample: Dict[str, Dict[str, Any]]) -> None:
        """Add five Trim Galore columns to the general statistics table."""

        gen_stats: Dict[str, Dict[str, float]] = {}
        for s_name, payload in data_by_sample.items():
            rp = payload.get("read_processing", {}) or {}
            bp = payload.get("basepair_processing", {}) or {}

            total_reads = rp.get("total_reads", 0) or 0
            reads_with_adapter = rp.get("reads_with_adapter", 0) or 0
            reads_written = rp.get("reads_written", 0) or 0
            total_bp = bp.get("total_bp_processed", 0) or 0
            quality_trimmed = bp.get("quality_trimmed_bp", 0) or 0

            row: Dict[str, float] = {
                "tg_total_reads": total_reads,
                "tg_pct_with_adapter": _safe_pct(reads_with_adapter, total_reads),
                "tg_pct_passing": _safe_pct(reads_written, total_reads),
                "tg_pct_quality_trimmed": _safe_pct(quality_trimmed, total_bp),
                "tg_total_bp_written": bp.get("total_bp_written", 0) or 0,
            }
            gen_stats[s_name] = row

        headers = {
            "tg_pct_with_adapter": {
                "title": "% Adapter",
                "description": "% reads where at least one adapter was detected",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "OrRd",
                "format": "{:,.1f}",
            },
            "tg_pct_passing": {
                "title": "% Pass",
                "description": "% reads passing all filters into trimmed output",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.1f}",
            },
            "tg_pct_quality_trimmed": {
                "title": "% Q-trim",
                "description": "% bases removed by 3'-end quality trimming",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "Blues",
                "format": "{:,.2f}",
            },
            "tg_total_reads": {
                "title": "Reads",
                "description": "Total reads processed",
                "scale": "Greys",
                "shared_key": "read_count",
                "hidden": True,
            },
            "tg_total_bp_written": {
                "title": "BP written",
                "description": "Total basepairs written to trimmed output",
                "scale": "BuGn",
                "shared_key": "base_count",
                "hidden": True,
            },
        }
        self.general_stats_addcols(gen_stats, headers)

    # ── Section: filtering disposition bar plot ─────────────────────────

    def _filtered_reads_plot(self, data_by_sample: Dict[str, Dict[str, Any]]):
        """Stacked bar plot of read disposition: passing / too_short /
        too_long / too_many_n / discarded_untrimmed."""
        bar_data: Dict[str, Dict[str, int]] = {}
        for s_name, payload in data_by_sample.items():
            rp = payload.get("read_processing", {}) or {}
            bar_data[s_name] = {
                "passing": rp.get("reads_written", 0) or 0,
                "too_short": rp.get("reads_too_short", 0) or 0,
                "too_long": rp.get("reads_too_long", 0) or 0,
                "too_many_n": rp.get("reads_too_many_n", 0) or 0,
                "discarded_untrimmed": rp.get("reads_discarded_untrimmed", 0) or 0,
            }
        cats = {
            "passing": {"name": "Passed filters", "color": "#7cb5ec"},
            "too_short": {"name": "Too short", "color": "#f7a35c"},
            "too_long": {"name": "Too long", "color": "#90ed7d"},
            "too_many_n": {"name": "Too many Ns", "color": "#8085e9"},
            "discarded_untrimmed": {
                "name": "Discarded (untrimmed)",
                "color": "#e4d354",
            },
        }
        return bargraph.plot(
            bar_data,
            cats,
            {
                "id": "trim_galore_filtered_reads_plot",
                "title": "Trim Galore: Read Filtering",
                "ylab": "Number of reads",
                "cpswitch_counts_label": "Number of reads",
            },
        )

    # ── Section: adapter length distribution ────────────────────────────

    def _adapter_length_plot(self, data_by_sample: Dict[str, Dict[str, Any]]):
        """Per-sample adapter length distribution. When a sample has multiple
        adapters (rare), each adapter becomes its own trace under a
        composite key."""
        line_data: Dict[str, Dict[int, int]] = {}
        for s_name, payload in data_by_sample.items():
            adapters = payload.get("adapter_trimming") or []
            if len(adapters) == 1:
                a = adapters[0]
                dist = a.get("length_distribution") or {}
                line_data[s_name] = {int(k): v for k, v in dist.items()}
            else:
                for idx, a in enumerate(adapters, start=1):
                    name = a.get("name") or f"adapter_{idx}"
                    composite = f"{s_name}: {name}"
                    dist = a.get("length_distribution") or {}
                    line_data[composite] = {int(k): v for k, v in dist.items()}
        if not line_data:
            return None
        return linegraph.plot(
            line_data,
            {
                "id": "trim_galore_adapter_length_plot",
                "title": "Trim Galore: Adapter Length Distribution",
                "xlab": "Adapter overlap length (bp)",
                "ylab": "Reads (count)",
                "ymin": 0,
                "tt_label": "{point.x} bp: {point.y:,} reads",
            },
        )


# ── Helpers ─────────────────────────────────────────────────────────────────


_FASTQ_SUFFIX_RE = re.compile(r"\.(fastq|fq)(\.gz)?$", re.IGNORECASE)


def _strip_fastq_suffix(name: str) -> str:
    """Strip `.fastq[.gz]` / `.fq[.gz]` from the input filename so the sample
    name is the bare basename (`sample_R1` rather than `sample_R1.fastq.gz`)."""
    return _FASTQ_SUFFIX_RE.sub("", name)


def _safe_pct(part: int, total: int) -> float:
    if not total:
        return 0.0
    return 100.0 * part / total
