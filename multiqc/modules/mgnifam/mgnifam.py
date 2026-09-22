import json
import logging
import math
from collections import Counter
from typing import Any

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph
from multiqc.plots.table_object import ColumnDict, ValueT
from multiqc.types import ColumnKey, SampleName, SectionAlert

log = logging.getLogger(__name__)

SCHEMA_VERSION = 1
COMMANDS = ("generate_families", "update_families")
STATS_SUFFIX = "_stats.json"
# Retention is stored exactly; for display it is binned into twentieths (0.05 wide).
RETENTION_BINS = 20


class MultiqcModule(BaseMultiqcModule):
    """
    mgnifam builds protein families from clusters of sequences by iterative HMM search
    (`mgnifam generate_families`), and refreshes existing family models against a newer
    sequence database (`mgnifam update_families`).

    The module parses the run summary each command writes after a completed run:
    `<chunk>_stats.json` for `generate_families` and `<chunk>_updated_stats.json` for
    `update_families`. One file describes one chunk, and each chunk is one sample. The
    sample name is the file name without the `_stats.json` suffix, so a chunk `1` reports
    as `1` for generation and `1_updated` for an update. Files with a `schema_version`
    other than 1 are skipped with a warning.

    A summary is only written when a run completes, so every sample in the report comes
    from consumable output. A run that completed with contained family crashes (exit
    status 3) is flagged in the outcomes section.

    Version 3.1.0 of mgnifam is tested.
    """

    def __init__(self):
        super().__init__(
            name="mgnifam",
            anchor="mgnifam",
            href="https://github.com/vagkaratzas/mgnifam",
            info="Iterative HMM-based protein family generation over very large sequence databases.",
            doi="10.5281/zenodo.22879938",
            license="MIT License",
            license_url="https://github.com/vagkaratzas/mgnifam/blob/main/LICENSE",
        )

        self.mgnifam_data: dict[str, dict[str, Any]] = {}
        for f in self.find_log_files("mgnifam"):
            try:
                stats = json.loads(f["f"])
            except ValueError as error:
                log.warning(f"Skipping {f['fn']}: not valid JSON ({error})")
                continue
            # The version gate decides whether the documented fields can be trusted at all.
            schema_version = stats.get("schema_version") if isinstance(stats, dict) else None
            if schema_version != SCHEMA_VERSION:
                log.warning(f"Skipping {f['fn']}: unsupported mgnifam schema_version {schema_version}")
                continue
            if stats["command"] not in COMMANDS:
                log.warning(f"Skipping {f['fn']}: unknown mgnifam command {stats['command']}")
                continue
            s_name = self.clean_s_name(f["fn"].removesuffix(STATS_SUFFIX), f)
            if s_name in self.mgnifam_data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, s_name=s_name)
            self.add_software_version(stats["version"], s_name)
            self.mgnifam_data[s_name] = stats

        self.mgnifam_data = self.ignore_samples(self.mgnifam_data)
        if len(self.mgnifam_data) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.mgnifam_data)} reports")

        self._add_general_stats()
        self._add_outcomes_section()
        self._add_histogram_section(
            "full_msa_size",
            "Full MSA size",
            "Sequences in the full alignment",
            "Each successful family's full alignment holds every sequence its final model recruited. "
            "Small alignments mark narrow families; very large ones can mark a family that absorbed "
            "unrelated sequences.",
        )
        self._add_histogram_section(
            "model_length",
            "Model length",
            "Match states in the family model",
            "The length of each successful family's final HMM, in match states.",
        )
        self._add_histogram_section(
            "model_length_change",
            "Model length change",
            "Match states gained or lost by the update",
            "Final model length minus the input model length, for every family that reached a final "
            "model. With "
            "`--skip_refine` models are not rebuilt, so every family sits at 0. In refine mode a model "
            "can shrink as columns are trimmed, or grow when new recruits support extra columns.",
            update_only=True,
        )
        self._add_histogram_section(
            "retention",
            "Retention",
            "Fraction of round 1 recruits still present at the end",
            "The share of the sequences recruited in the first search round that are still members "
            "when the family finishes. Under `--skip_refine` it is 1.0 by construction. In refine mode, "
            "values near 1 mean the family stayed stable. Values "
            "are grouped into bins 0.05 wide for display; `multiqc_mgnifam_histograms.json` keeps them exact.",
            update_only=True,
        )

        self.write_data_file(
            {
                s_name: {
                    "command": stats["command"],
                    "chunk_id": stats["chunk_id"],
                    "version": stats["version"],
                    "exit_status": stats["exit_status"],
                    **stats["families"],
                }
                for s_name, stats in self.mgnifam_data.items()
            },
            "multiqc_mgnifam",
        )
        # Exact values, including the histograms the report does not plot.
        self.write_data_file(
            {s_name: stats["histograms"] for s_name, stats in self.mgnifam_data.items()},
            "multiqc_mgnifam_histograms",
            data_format="json",
        )

    def _add_general_stats(self):
        data: dict[SampleName | str, dict[ColumnKey | str, ValueT]] = {}
        for s_name, stats in self.mgnifam_data.items():
            families = stats["families"]
            row: dict[ColumnKey | str, ValueT] = {
                "families_input": families["input"],
                "converged": families["converged"],
            }
            # An empty chunk is valid output: leave the percentage out rather than divide by zero.
            if families["input"]:
                row["successful_pct"] = 100 * families["successful"] / families["input"]
            retention = stats["histograms"].get("retention", {})
            if retention:
                row["mean_retention"] = sum(float(value) * count for value, count in retention.items()) / sum(
                    retention.values()
                )
            data[s_name] = row

        headers: dict[str, ColumnDict] = {
            "families_input": {
                "title": "Families In",
                "description": "Clusters (generation) or models (update) processed",
                "scale": "Blues",
                "format": "{:,.0f}",
            },
            "successful_pct": {
                "title": "% Successful",
                "description": "Percentage of input families that produced a family",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "converged": {
                "title": "Converged",
                "description": "Successful families that converged before the round limit",
                "scale": "Greens",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "mean_retention": {
                "title": "Retention",
                "description": "Mean fraction of round 1 recruits still present at the end of an update",
                "min": 0,
                "max": 1,
                "scale": "YlGnBu",
                "format": "{:,.2f}",
            },
        }
        general_stats_headers = self.get_general_stats_headers(all_headers=headers)
        if general_stats_headers:
            self.general_stats_addcols(data, general_stats_headers)

    def _add_outcomes_section(self):
        data = {
            s_name: {"Successful": stats["families"]["successful"], **stats["discard_reasons"]}
            for s_name, stats in self.mgnifam_data.items()
        }
        categories = ["Successful"] + sorted(
            {reason for stats in self.mgnifam_data.values() for reason in stats["discard_reasons"]}
        )
        crashed = sorted(s_name for s_name, stats in self.mgnifam_data.items() if stats["families"]["crashed"])
        self.add_section(
            name="Family outcomes",
            anchor="mgnifam-outcomes",
            description="Families that were generated or updated successfully, and the reasons the rest were discarded.",
            helptext="""
Each input cluster (for `generate_families`) or model (for `update_families`) ends in exactly
one outcome: it becomes a successful family, or it is discarded with the reason shown here.

A family that failed on an internal error is contained and recorded as a discard, so the chunk
still completes, but mgnifam then exits with status 3. Chunks where that happened are listed
in the warning above the plot.
""",
            plot=bargraph.plot(
                data,
                {category: {"name": category.capitalize()} for category in categories},
                pconfig={
                    "id": "mgnifam-outcomes-plot",
                    "title": "mgnifam: Family outcomes",
                    "ylab": "Families",
                },
            ),
            alerts=SectionAlert(
                message=(
                    f"**{len(crashed)} chunk{'s' if len(crashed) != 1 else ''}** completed with "
                    "families that failed on an internal error and were recorded as discards."
                ),
                level="warning",
                affected_samples=crashed,
            )
            if crashed
            else None,
        )

    def _add_histogram_section(self, key: str, name: str, xlab: str, helptext: str, update_only: bool = False):
        samples = {
            s_name: stats
            for s_name, stats in self.mgnifam_data.items()
            if not update_only or stats["command"] == "update_families"
        }
        if not samples:
            return
        data: dict[str, dict[float, int]] = {}
        for s_name, stats in samples.items():
            counts: Counter[float] = Counter()
            for value, count in stats["histograms"][key].items():
                x = float(value)
                if key == "retention":
                    x = math.floor(x * RETENTION_BINS) / RETENTION_BINS
                counts[x] += count
            if counts:
                data[s_name] = dict(sorted(counts.items()))
        empty = sorted(set(samples) - set(data))
        anchor = f"mgnifam-{key.replace('_', '-')}"
        self.add_section(
            name=name,
            anchor=anchor,
            description=f"Number of families per value of: {xlab.lower()}.",
            helptext=helptext,
            plot=linegraph.plot(
                data,
                pconfig={
                    "id": f"{anchor}-plot",
                    "title": f"mgnifam: {name}",
                    "xlab": xlab,
                    "ylab": "Families",
                    "y_decimals": False,
                },
            )
            if data
            else None,
            alerts=SectionAlert(
                message=(
                    f"**{len(empty)} chunk{'s' if len(empty) != 1 else ''}** "
                    "with no families to plot hidden from this plot."
                ),
                level="info",
                affected_samples=empty,
            )
            if empty
            else None,
        )
