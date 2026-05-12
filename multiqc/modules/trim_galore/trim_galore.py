import json
import logging
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple, cast

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound, SampleGroupingConfig
from multiqc.plots import bargraph, linegraph, table
from multiqc.types import ColumnKey, SampleName, SectionAlert

log = logging.getLogger(__name__)


SCHEMA_MAJOR_SUPPORTED = 1


class MultiqcModule(BaseMultiqcModule):
    """
    [Trim Galore](https://github.com/FelixKrueger/TrimGalore) provides
    consistent quality and adapter trimming for next-generation sequencing
    data, with special handling for Reduced Representation Bisulfite
    Sequencing (RRBS) and small-RNA libraries.

    This MultiQC module supports Trim Galore v2.0, which is a Rust
    rewrite of the original Perl-based v0.6 that has a new JSON output
    file summarising results. The earlier v0.6 versions of Trim Galore
    that wrapped Cutadapt are supported with reporting via the Cutadapt module.
    The old log format is still produced, but the Cutadapt module search pattern
    is configured to skip reports mentioning Trim Galore v2+. If you delete the
    JSON but keep the v2 text file, the sample will not be reported by either module.

    #### Paired-end sample grouping

    R1 and R2 of a paired-end sample are grouped automatically into a single
    row in the General Statistics and Pair Validation tables. Click the
    expand arrow on a grouped row to see the per-read values. The grouping
    is derived from the file list inside each JSON report, so it does not
    depend on filename patterns.

    To disable automatic grouping and show one row per read everywhere:

    ```yaml
    trim_galore_config:
      auto_group_pairs: false
    ```

    Auto-grouping can be combined with the global
    [`table_sample_merge`](../reports/customisation.md#sample-grouping)
    config option to merge further — for example to group lanes of an
    already-paired sample.
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
            doi="10.5281/zenodo.5127898",
        )

        data_by_sample: Dict[str, Dict[str, Any]] = {}
        self._pair_key_to_samples: Dict[Tuple[str, ...], List[str]] = {}
        for f in self.find_log_files("trim_galore", filehandles=True):
            parsed = self._parse_log(f)
            if parsed is None:
                continue
            s_name, payload, pair_key = parsed
            if self.is_ignore_sample(s_name):
                continue
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            data_by_sample[s_name] = payload
            self._pair_key_to_samples.setdefault(pair_key, []).append(s_name)
            self.add_data_source(f, s_name=s_name)
            self.add_software_version(payload["trim_galore_version"], s_name)

        if not data_by_sample:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(data_by_sample)} reports")

        auto_group_pairs = getattr(config, "trim_galore_config", {}).get("auto_group_pairs", True)
        self._pair_display_by_key: Dict[Tuple[str, ...], str] = {}
        explicit_groups: Optional[Dict[str, List[str]]] = None
        if auto_group_pairs:
            explicit_groups, self._pair_display_by_key = self._derive_auto_groups(data_by_sample)

        self._general_stats_table(data_by_sample, explicit_groups)

        self.add_section(
            name="Filtered Reads",
            anchor="trim_galore_filtered_reads",
            description="Number of reads after trimming.",
            helptext=(
                "The `passing` bucket reflects reads that survived quality trimming, adapter trimming and N-content filtering "
                "(i.e. what made it into the trimmed output); the remaining buckets show which filter rejected the rest."
            ),
            plot=self._filtered_reads_plot(data_by_sample),
        )

        adapter_plot = self._adapter_length_plot(data_by_sample)
        if adapter_plot is not None:
            self.add_section(
                name="Adapter Length Distribution",
                anchor="trim_galore_adapter_lengths",
                description="Histogram of adapter match lengths per read.",
                helptext=(
                    "Note that a tall left tail at 1 bp is normal: Trim Galore's default `--stringency 1` accepts "
                    "single-base overlaps, which represent random hits rather than real adapter contamination."
                ),
                plot=adapter_plot,
            )

        pair_validation_plot, pv_dropped = self._pair_validation_plot(data_by_sample)
        if pair_validation_plot is not None:
            self.add_section(
                name="Pair Validation",
                anchor="trim_galore_pair_validation",
                description="Outcomes of paired-end validation",
                helptext=(
                    "Shows pairs analysed, pairs removed (broken down by reason), and reads left unpaired after a partner "
                    "was discarded. R1 and R2 of a pair are collapsed into a single row (the data is identical between them)."
                ),
                plot=pair_validation_plot,
                alerts=SectionAlert(
                    message=(
                        f"**{len(pv_dropped)} sample{'s' if len(pv_dropped) != 1 else ''}** with less than 0.1% "
                        "of pairs affected hidden from this table."
                    ),
                    affected_samples=pv_dropped,
                )
                if pv_dropped
                else None,
            )

        poly_plot, poly_dropped = self._poly_trimming_plot(data_by_sample)
        if poly_plot is not None:
            self.add_section(
                name="Poly-A / Poly-G Trimming",
                anchor="trim_galore_poly_trimming",
                description="Reads and bases removed by Trim Galore's poly-A and poly-G tail trimmers.",
                plot=poly_plot,
                alerts=SectionAlert(
                    message=(
                        f"**{len(poly_dropped)} sample{'s' if len(poly_dropped) != 1 else ''}** with no "
                        "poly-A/G trimming hidden from this table."
                    ),
                    affected_samples=poly_dropped,
                )
                if poly_dropped
                else None,
            )

        rrbs_plot, rrbs_dropped = self._rrbs_plot(data_by_sample)
        if rrbs_plot is not None:
            self.add_section(
                name="RRBS Trimming",
                anchor="trim_galore_rrbs",
                description="Reduced Representation Bisulfite Sequencing (RRBS) end-repair trimming counts.",
                plot=rrbs_plot,
                alerts=SectionAlert(
                    message=(
                        f"**{len(rrbs_dropped)} sample{'s' if len(rrbs_dropped) != 1 else ''}** with no "
                        "RRBS trimming hidden from this table."
                    ),
                    affected_samples=rrbs_dropped,
                )
                if rrbs_dropped
                else None,
            )

        self.write_data_file(_flatten_for_data_file(data_by_sample), "multiqc_trim_galore")

    def _parse_log(self, f) -> Optional[Tuple[str, Dict[str, Any], Tuple[str, ...]]]:
        try:
            payload = json.load(f["f"])
        except json.JSONDecodeError as e:
            log.warning(f"Could not parse {f['fn']!r}: {e}")
            return None

        # `tool` check uses .get() because we don't yet know this is a Trim
        # Galore JSON — could be a foreign JSON that matched the filename glob.
        if payload.get("tool") != "Trim Galore":
            log.debug(f"Skipping {f['fn']!r}: tool field is {payload.get('tool')!r}, not Trim Galore")
            return None

        schema = payload["schema_version"]
        major = int(str(schema).split(".")[0])
        if major != SCHEMA_MAJOR_SUPPORTED:
            log.error(
                f"Skipping {f['fn']!r}: schema_version {schema!r} has incompatible major version "
                f"(expected major v{SCHEMA_MAJOR_SUPPORTED}). This module needs updating to handle the new schema."
            )
            return None

        # PE reports list both filenames in `input_filenames`; pick the one for
        # this read_number so R1 and R2 produce distinct samples. SE reports
        # have read_number=1 and a 1-element list, so the same logic works.
        # `input_filenames` is byte-identical across R1 and R2 of the same pair,
        # so it doubles as a stable pair key.
        input_filenames = payload["input_filenames"]
        read_number = payload["read_number"]
        s_name = self.clean_s_name(input_filenames[read_number - 1], f)
        return s_name, payload, tuple(input_filenames)

    def _derive_auto_groups(
        self,
        data_by_sample: Dict[str, Dict[str, Any]],
    ) -> Tuple[Dict[str, List[str]], Dict[Tuple[str, ...], str]]:
        """Build module-supplied groups from tool-derived pair info.

        Returns `(explicit_groups, pair_display_by_key)`. When the user has
        `config.table_sample_merge` set, name patterns layer on top: each
        auto-group's display name is run through `groups_for_sample` to
        merge into a super-group.
        """
        auto_groups: Dict[str, List[str]] = {}
        pair_display_by_key: Dict[Tuple[str, ...], str] = {}
        for pair_key, members in self._pair_key_to_samples.items():
            sorted_members = sorted(members)
            # _clean_fastq_pair needs two names; SE pairs have only one.
            if len(sorted_members) >= 2:
                display = self._clean_fastq_pair(sorted_members[0], sorted_members[1]) or "_".join(sorted_members)
            else:
                display = sorted_members[0]
            auto_groups[display] = sorted_members
            pair_display_by_key[pair_key] = display

        if not config.table_sample_merge:
            return auto_groups, pair_display_by_key

        by_super: Dict[str, List[str]] = defaultdict(list)
        for auto_display, samples in auto_groups.items():
            super_name, _ = self.groups_for_sample(SampleName(auto_display))
            by_super[str(super_name)].extend(samples)
        return dict(by_super), pair_display_by_key

    def _general_stats_table(
        self,
        data_by_sample: Dict[str, Dict[str, Any]],
        explicit_groups: Optional[Dict[str, List[str]]] = None,
    ) -> None:
        gen_stats: Dict[str, Dict[ColumnKey, Any]] = {}
        for s_name, payload in data_by_sample.items():
            rp = payload["read_processing"]
            bp = payload["basepair_processing"]
            total = rp["total_reads"]
            total_bp = bp["total_bp_processed"]
            gen_stats[s_name] = {
                ColumnKey("tg_total_reads"): total,
                ColumnKey("tg_pct_with_adapter"): _safe_pct(rp["reads_with_adapter"], total),
                ColumnKey("tg_pct_passing"): _safe_pct(rp["reads_written"], total),
                ColumnKey("tg_pct_quality_trimmed"): _safe_pct(bp["quality_trimmed_bp"], total_bp),
                ColumnKey("tg_total_bp_written"): bp["total_bp_written"],
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
            cast(Any, gen_stats),
            headers,
            group_samples_config=SampleGroupingConfig(
                explicit_groups=explicit_groups,
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
            rp = payload["read_processing"]
            bar_data[s_name] = {
                "passing": rp["reads_written"],
                "too_short": rp["reads_too_short"],
                "too_long": rp["reads_too_long"],
                "too_many_n": rp["reads_too_many_n"],
                "discarded_untrimmed": rp["reads_discarded_untrimmed"],
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
            adapters = payload["adapter_trimming"]
            for a in adapters:
                key = f"{s_name} ({a['name']})" if len(adapters) > 1 else s_name
                line_data[key] = {int(k): v for k, v in a["length_distribution"].items()}
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

    def _pair_validation_plot(self, data_by_sample: Dict[str, Dict[str, Any]]):
        # pair_validation is identical between R1 and R2 JSONs of the same pair —
        # collapse them by the tool-derived pair_key. Skip SE samples (no pv).
        all_rows: Dict[str, Dict[str, int]] = {}
        for pair_key, members in self._pair_key_to_samples.items():
            pv = data_by_sample[members[0]]["pair_validation"]
            if not pv:  # SE samples — `pair_validation` is null
                continue
            display = self._pair_display_by_key.get(pair_key, members[0])
            all_rows[display] = {
                "pairs_analyzed": pv["pairs_analyzed"],
                "pairs_removed": pv["pairs_removed"],
                "pairs_removed_n": pv["pairs_removed_n"],
                "pairs_removed_too_long": pv["pairs_removed_too_long"],
                "r1_unpaired": pv["r1_unpaired"],
                "r2_unpaired": pv["r2_unpaired"],
            }

        # `pairs_removed_*` are sub-reasons of `pairs_removed` so are not added
        # to the affected-fraction sum here (would double-count).
        def _kept(r: Dict[str, int]) -> bool:
            return bool(r["pairs_analyzed"]) and (
                (r["pairs_removed"] + r["r1_unpaired"] + r["r2_unpaired"]) / r["pairs_analyzed"] > 0.001
            )

        rows = {n: r for n, r in all_rows.items() if _kept(r)}
        dropped = sorted(n for n, r in all_rows.items() if not _kept(r))
        if not rows:
            return None, dropped
        headers: Dict[str, Dict[str, Any]] = {
            "pairs_analyzed": {
                "title": "Pairs analysed",
                "description": "Total read pairs examined by pair validation",
                "scale": "Greys",
                "shared_key": "read_count",
            },
            "pairs_removed": {
                "title": "Pairs removed",
                "description": "Total read pairs removed by pair validation",
                "scale": "OrRd",
                "shared_key": "read_count",
            },
            "pairs_removed_n": {
                "title": "Pairs removed (N content)",
                "description": "Read pairs removed because they exceeded the N-content threshold",
                "scale": "OrRd",
                "shared_key": "read_count",
            },
            "pairs_removed_too_long": {
                "title": "Pairs removed (too long)",
                "description": "Read pairs removed because one or both reads exceeded the maximum length",
                "scale": "OrRd",
                "shared_key": "read_count",
            },
            "r1_unpaired": {
                "title": "R1 unpaired",
                "description": "R1 reads left unpaired after their R2 partner was discarded",
                "scale": "Oranges",
                "shared_key": "read_count",
            },
            "r2_unpaired": {
                "title": "R2 unpaired",
                "description": "R2 reads left unpaired after their R1 partner was discarded",
                "scale": "Oranges",
                "shared_key": "read_count",
            },
        }
        return (
            table.plot(
                rows,
                cast(Any, headers),
                pconfig={
                    "namespace": self.name,
                    "id": "trim_galore_pair_validation_table",
                    "title": "Trim Galore: Pair Validation",
                },
            ),
            dropped,
        )

    def _poly_trimming_plot(self, data_by_sample: Dict[str, Dict[str, Any]]):
        rows: Dict[str, Dict[str, int]] = {}
        dropped: List[str] = []
        for s_name, payload in data_by_sample.items():
            pa = payload["poly_a_trimming"]
            pg = payload["poly_g_trimming"]
            row = {
                "poly_a_reads_trimmed": pa["reads_trimmed"],
                "poly_a_bases_removed": pa["bases_removed"],
                "poly_g_reads_trimmed": pg["reads_trimmed"],
                "poly_g_bases_removed": pg["bases_removed"],
            }
            if any(row.values()):
                rows[s_name] = row
            else:
                dropped.append(s_name)
        if not rows:
            return None, sorted(dropped)
        headers: Dict[str, Dict[str, Any]] = {
            "poly_a_reads_trimmed": {
                "title": "Poly-A reads trimmed",
                "description": "Reads with a poly-A tail trimmed",
                "scale": "Purples",
                "shared_key": "read_count",
            },
            "poly_a_bases_removed": {
                "title": "Poly-A bp removed",
                "description": "Bases removed by poly-A trimming",
                "scale": "Purples",
                "shared_key": "base_count",
            },
            "poly_g_reads_trimmed": {
                "title": "Poly-G reads trimmed",
                "description": "Reads with a poly-G tail trimmed",
                "scale": "Greens",
                "shared_key": "read_count",
            },
            "poly_g_bases_removed": {
                "title": "Poly-G bp removed",
                "description": "Bases removed by poly-G trimming",
                "scale": "Greens",
                "shared_key": "base_count",
            },
        }
        return (
            table.plot(
                rows,
                cast(Any, headers),
                pconfig={
                    "namespace": self.name,
                    "id": "trim_galore_poly_trimming_table",
                    "title": "Trim Galore: Poly-A / Poly-G Trimming",
                },
            ),
            sorted(dropped),
        )

    def _rrbs_plot(self, data_by_sample: Dict[str, Dict[str, Any]]):
        rows: Dict[str, Dict[str, int]] = {}
        dropped: List[str] = []
        for s_name, payload in data_by_sample.items():
            rr = payload["rrbs"]
            row = {
                "rrbs_trimmed_3prime": rr["trimmed_3prime"],
                "rrbs_trimmed_5prime": rr["trimmed_5prime"],
                "rrbs_r2_clipped_5prime": rr["r2_clipped_5prime"],
            }
            if any(row.values()):
                rows[s_name] = row
            else:
                dropped.append(s_name)
        if not rows:
            return None, sorted(dropped)
        headers: Dict[str, Dict[str, Any]] = {
            "rrbs_trimmed_3prime": {
                "title": "3' trimmed",
                "description": "Reads trimmed at the 3' end for RRBS end-repair",
                "scale": "Blues",
                "shared_key": "read_count",
            },
            "rrbs_trimmed_5prime": {
                "title": "5' trimmed",
                "description": "Reads trimmed at the 5' end for RRBS end-repair",
                "scale": "Blues",
                "shared_key": "read_count",
            },
            "rrbs_r2_clipped_5prime": {
                "title": "R2 5' clipped",
                "description": "R2 reads with the 5' end clipped (directional RRBS)",
                "scale": "Blues",
                "shared_key": "read_count",
            },
        }
        return (
            table.plot(
                rows,
                cast(Any, headers),
                pconfig={
                    "namespace": self.name,
                    "id": "trim_galore_rrbs_table",
                    "title": "Trim Galore: RRBS Trimming",
                },
            ),
            sorted(dropped),
        )


def _flatten_for_data_file(data_by_sample: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    flat: Dict[str, Dict[str, Any]] = {}
    for s_name, p in data_by_sample.items():
        rp = p["read_processing"]
        bp = p["basepair_processing"]
        pa = p["poly_a_trimming"]
        pg = p["poly_g_trimming"]
        rr = p["rrbs"]
        # `pair_validation` is explicitly null for SE samples — keep `or {}`
        # so the .get() calls below produce None for the SE pair columns.
        pv = p["pair_validation"] or {}
        flat[s_name] = {
            "trim_galore_version": p["trim_galore_version"],
            "mode": p["mode"],
            "total_reads": rp["total_reads"],
            "reads_with_adapter": rp["reads_with_adapter"],
            "reads_written": rp["reads_written"],
            "reads_too_short": rp["reads_too_short"],
            "reads_too_long": rp["reads_too_long"],
            "reads_too_many_n": rp["reads_too_many_n"],
            "reads_discarded_untrimmed": rp["reads_discarded_untrimmed"],
            "total_bp_processed": bp["total_bp_processed"],
            "total_bp_written": bp["total_bp_written"],
            "quality_trimmed_bp": bp["quality_trimmed_bp"],
            "pairs_analyzed": pv.get("pairs_analyzed"),
            "pairs_removed": pv.get("pairs_removed"),
            "pairs_removed_n": pv.get("pairs_removed_n"),
            "pairs_removed_too_long": pv.get("pairs_removed_too_long"),
            "r1_unpaired": pv.get("r1_unpaired"),
            "r2_unpaired": pv.get("r2_unpaired"),
            "poly_a_reads_trimmed": pa["reads_trimmed"],
            "poly_a_bases_removed": pa["bases_removed"],
            "poly_g_reads_trimmed": pg["reads_trimmed"],
            "poly_g_bases_removed": pg["bases_removed"],
            "rrbs_trimmed_3prime": rr["trimmed_3prime"],
            "rrbs_trimmed_5prime": rr["trimmed_5prime"],
            "rrbs_r2_clipped_5prime": rr["r2_clipped_5prime"],
        }
    return flat


def _safe_pct(part: int, total: int) -> float:
    if not total:
        return 0.0
    return 100.0 * part / total
