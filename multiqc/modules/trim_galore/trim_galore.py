import json
import logging
from typing import Any, Dict, Optional, Tuple, Union

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound, SampleGroupingConfig
from multiqc.plots import bargraph, linegraph
from multiqc.types import ColumnKey, SampleName

log = logging.getLogger(__name__)


SCHEMA_VERSION_SUPPORTED = 1
TOOL_NAME = "Trim Galore"


class MultiqcModule(BaseMultiqcModule):
    """
    Parses the structured `*_trimming_report.json` (schema v1) emitted by
    Trim Galore v2.x alongside the legacy text report.

    Trim Galore v2.x also writes a `*_trimming_report.txt` carrying a
    "This is cutadapt ..." backwards-compatibility shim. The `cutadapt`
    module's search pattern excludes v2.x text reports (matched on the
    `Trim Galore version:` header) so samples are not double-counted.
    Legacy Trim Galore v0.x / v1.x text reports continue to be parsed by
    the `cutadapt` module — they have no JSON sibling. If the JSON is
    removed but the v2.x text file is kept, the sample will not be
    reported by either module.
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
                "This module parses the structured JSON report (`*_trimming_report.json`) "
                "introduced in Trim Galore v2.0 alongside the legacy text report."
            ),
            doi="10.5281/zenodo.5127899",
        )

        data_by_sample: Dict[str, Dict[str, Any]] = {}
        for f in self.find_log_files("trim_galore", filehandles=True):
            parsed = self._parse_log(f)
            if parsed is None:
                continue
            s_name, payload = parsed
            if self.is_ignore_sample(s_name):
                continue
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            data_by_sample[s_name] = payload
            self.add_data_source(f, s_name=s_name)
            self.add_software_version(payload.get("trim_galore_version"), s_name)

        if not data_by_sample:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(data_by_sample)} reports")

        self._general_stats_table(data_by_sample)

        self.add_section(
            name="Filtered Reads",
            anchor="trim_galore_filtered_reads",
            description=(
                "Disposition of reads after trimming. The `passing` bucket reflects reads "
                "that survived quality trimming, adapter trimming and N-content filtering "
                "(i.e. what made it into the trimmed output); the remaining buckets show "
                "which filter rejected the rest."
            ),
            plot=self._filtered_reads_plot(data_by_sample),
        )

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

        self.write_data_file(_flatten_for_data_file(data_by_sample), "multiqc_trim_galore")

    def _parse_log(self, f) -> Optional[Tuple[str, Dict[str, Any]]]:
        try:
            payload = json.load(f["f"])
        except json.JSONDecodeError as e:
            log.warning(f"Could not parse {f['fn']!r}: {e}")
            return None

        if payload.get("tool") != TOOL_NAME:
            log.debug(f"Skipping {f['fn']!r}: tool field is {payload.get('tool')!r}, not {TOOL_NAME!r}")
            return None
        schema = payload.get("schema_version")
        if schema != SCHEMA_VERSION_SUPPORTED:
            log.error(
                f"Skipping {f['fn']!r}: unsupported schema_version {schema!r} (expected {SCHEMA_VERSION_SUPPORTED}). "
                "This module needs updating to handle the new schema."
            )
            return None

        # PE reports list BOTH input filenames in `input_filenames`; use
        # `read_number` to pick the matching one so R1 and R2 produce
        # distinct samples. SE reports have `read_number: 1` and a
        # single-entry list, so the same logic gives the correct SE behaviour.
        input_filenames = payload.get("input_filenames") or []
        read_number = payload.get("read_number") or 1
        if not input_filenames:
            log.warning(f"Skipping {f['fn']!r}: no input_filenames")
            return None
        idx = max(0, min(read_number - 1, len(input_filenames) - 1))
        s_name = self.clean_s_name(input_filenames[idx], f)
        return s_name, payload

    def _general_stats_table(self, data_by_sample: Dict[str, Dict[str, Any]]) -> None:
        gen_stats: Dict[Union[SampleName, str], Dict[Union[ColumnKey, str], Union[int, float, str, bool]]] = {}
        for s_name, payload in data_by_sample.items():
            rp = payload.get("read_processing", {}) or {}
            bp = payload.get("basepair_processing", {}) or {}
            total = rp.get("total_reads", 0) or 0
            total_bp = bp.get("total_bp_processed", 0) or 0
            gen_stats[s_name] = {
                ColumnKey("tg_total_reads"): total,
                ColumnKey("tg_pct_with_adapter"): _safe_pct(rp.get("reads_with_adapter", 0) or 0, total),
                ColumnKey("tg_pct_passing"): _safe_pct(rp.get("reads_written", 0) or 0, total),
                ColumnKey("tg_pct_quality_trimmed"): _safe_pct(bp.get("quality_trimmed_bp", 0) or 0, total_bp),
                ColumnKey("tg_total_bp_written"): bp.get("total_bp_written", 0) or 0,
            }

        headers: Dict[str, Dict[str, Any]] = {
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
                "format": "{:,.1f}",
            },
            "tg_total_reads": {
                "title": "Reads",
                "description": "Total reads processed",
                "scale": "Greys",
                "shared_key": "read_count",
            },
            "tg_total_bp_written": {
                "title": "BP written",
                "description": "Total basepairs written to trimmed output",
                "scale": "BuGn",
                "shared_key": "base_count",
                "hidden": True,
            },
        }
        self.general_stats_addcols(
            gen_stats,
            headers,
            group_samples_config=SampleGroupingConfig(
                cols_to_sum=[
                    ColumnKey("tg_total_reads"),
                    ColumnKey("tg_total_bp_written"),
                ],
                cols_to_weighted_average=[
                    (ColumnKey("tg_pct_with_adapter"), ColumnKey("tg_total_reads")),
                    (ColumnKey("tg_pct_passing"), ColumnKey("tg_total_reads")),
                    (ColumnKey("tg_pct_quality_trimmed"), ColumnKey("tg_total_bp_written")),
                ],
            ),
        )

    def _filtered_reads_plot(self, data_by_sample: Dict[str, Dict[str, Any]]):
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
            "passing": {"name": "Passed filters"},
            "too_short": {"name": "Too short"},
            "too_long": {"name": "Too long"},
            "too_many_n": {"name": "Too many Ns"},
            "discarded_untrimmed": {"name": "Discarded (untrimmed)"},
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

    def _adapter_length_plot(self, data_by_sample: Dict[str, Dict[str, Any]]):
        line_data: Dict[str, Dict[int, int]] = {}
        for s_name, payload in data_by_sample.items():
            adapters = payload.get("adapter_trimming") or []
            if not adapters:
                continue
            for idx, a in enumerate(adapters, start=1):
                name = a.get("name") or f"adapter_{idx}"
                key = f"{s_name} ({name})" if len(adapters) > 1 else s_name
                dist = a.get("length_distribution") or {}
                line_data[key] = {int(k): v for k, v in dist.items()}
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


def _flatten_for_data_file(data_by_sample: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    flat: Dict[str, Dict[str, Any]] = {}
    for s_name, p in data_by_sample.items():
        rp = p.get("read_processing", {}) or {}
        bp = p.get("basepair_processing", {}) or {}
        pv = p.get("pair_validation") or {}
        flat[s_name] = {
            "trim_galore_version": p.get("trim_galore_version"),
            "mode": p.get("mode"),
            "total_reads": rp.get("total_reads"),
            "reads_with_adapter": rp.get("reads_with_adapter"),
            "reads_written": rp.get("reads_written"),
            "reads_too_short": rp.get("reads_too_short"),
            "reads_too_long": rp.get("reads_too_long"),
            "reads_too_many_n": rp.get("reads_too_many_n"),
            "reads_discarded_untrimmed": rp.get("reads_discarded_untrimmed"),
            "total_bp_processed": bp.get("total_bp_processed"),
            "total_bp_written": bp.get("total_bp_written"),
            "quality_trimmed_bp": bp.get("quality_trimmed_bp"),
            "pairs_analyzed": pv.get("pairs_analyzed") if pv else None,
            "pairs_removed": pv.get("pairs_removed") if pv else None,
        }
    return flat


def _safe_pct(part: int, total: int) -> float:
    if not total:
        return 0.0
    return 100.0 * part / total
