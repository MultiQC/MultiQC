import json
import logging
import os

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph

log = logging.getLogger(__name__)

# Summaries this module reads (pmdshodh 0.1.0 and later).
FORMAT_VERSIONS = (1, 2, 3, 4)
BASES = ("A", "C", "G", "T")


def percent(num: float, den: float) -> float | None:
    """`100 num / den`, or None when nothing was counted."""
    return 100.0 * num / den if den > 0 else None


def damage_curves(summary: dict) -> dict[str, dict[int, float]]:
    """
    Damage rates by position from each molecule end, in percent, as `pmdshodh report` computes
    them: 5' C>T, and 3' G>A for double-stranded or 3' C>T for single-stranded libraries. Each
    rate is the read base's count over all read bases at that reference base and position.
    Positions are numbered from 1 (the terminal base).
    """
    damage = summary["damage"]
    ss = summary["mode"]["library"] == "ss"
    shown = {"five_prime": ("C", "T"), "three_prime": ("C", "T") if ss else ("G", "A")}
    curves: dict[str, dict[int, float]] = {}
    for end, (ref, base) in shown.items():
        table = damage[end][ref]
        rates = {}
        for i in range(damage["range"]):
            rate = percent(table[base][i], sum(table[b][i] for b in BASES))
            if rate is not None:
                rates[i + 1] = rate
        curves[end] = rates
    return curves


def parse_summary(summary: dict) -> dict | None:
    """The values this module reports from one run summary, or None if it is not one it reads."""
    if summary.get("tool", {}).get("name") != "pmdshodh":
        return None
    if summary.get("format_version") not in FORMAT_VERSIONS:
        return None
    reads = summary["reads"]
    filtering = summary["subcommand"] == "filter"
    first = summary["damage"]["first_position"]
    parsed = {
        "version": summary["tool"]["version"],
        "subcommand": summary["subcommand"],
        "input": summary["input"],
        "library": summary["mode"]["library"],
        "treatment": summary["mode"]["treatment"],
        "pmdtools_compat": summary["mode"]["pmdtools_compat"],
        "threshold": summary["filters"].get("threshold"),
        "input_reads": reads["input"],
        "scored": reads["scored"],
        "kept": reads["kept"] if filtering else None,
        "removed_by_threshold": reads["removed_by_threshold"] if filtering else None,
        "no_informative_sites": reads["no_informative_sites"],
        "skipped": dict(reads["skipped"]),
        "pct_kept": percent(reads["kept"], reads["scored"]) if filtering else None,
        "first_c_to_t": 100.0 * first["c_to_t"] if first["n"] > 0 else None,
        "first_c_to_t_se": 100.0 * first["se"] if first["n"] > 0 else None,
        "first_c_to_t_n": first["n"],
    }
    hist = summary["score_histogram"]
    scored = reads["scored"]
    parsed["score_histogram"] = (
        {
            round(hist["min"] + (i + 0.5) * hist["bin_width"], 6): 100.0 * c / scored
            for i, c in enumerate(hist["counts"])
        }
        if scored > 0
        else {}
    )
    lengths = summary["length_histogram"]
    parsed["lengths_scored"] = {int(k): v for k, v in lengths["scored"].items()}
    parsed["lengths_kept"] = {int(k): v for k, v in (lengths["kept"] or {}).items()}
    parsed["damage"] = damage_curves(summary)
    return parsed


class MultiqcModule(BaseMultiqcModule):
    """
    The module reads the JSON run summaries that pmdshodh writes, with `damage -o`,
    `filter --summary` or `score --summary`. It finds them by their content, whatever
    their file name, and reads `format_version` 1 to 4 (pmdshodh 0.1.0 and later).

    ```bash
    pmdshodh filter --threshold 3 --summary sample.pmdshodh.json sample.bam -o sample.pmd3.bam
    ```

    Samples are named after the alignment file that pmdshodh read (the summary's `input`
    field), cleaned as MultiQC cleans file names.

    Damage rates are computed from the counts in the summary, as `pmdshodh report` does:
    5′ C→T, and 3′ G→A for double-stranded or 3′ C→T for single-stranded libraries.
    """

    def __init__(self):
        super().__init__(
            name="pmdshodh",
            anchor="pmdshodh",
            href="https://github.com/dpatel511/pmdshodh",
            info="Scores and filters ancient DNA reads by post-mortem damage (PMD score).",
            doi="10.1073/pnas.1318934111",
            license="MIT License",
            license_url="https://github.com/dpatel511/pmdshodh/blob/main/LICENSE",
        )

        data: dict[str, dict] = {}
        for f in self.find_log_files("pmdshodh"):
            try:
                summary = json.loads(f["f"])
            except json.JSONDecodeError as e:
                log.warning(f"Could not parse {f['fn']} as JSON: {e}")
                continue
            try:
                parsed = parse_summary(summary)
            except (KeyError, TypeError) as e:
                log.warning(f"Skipping {f['fn']}: missing or invalid field {e}")
                continue
            if parsed is None:
                log.debug(f"Skipping {f['fn']}: not a pmdshodh run summary this module reads")
                continue
            s_name = self.clean_s_name(os.path.basename(parsed["input"]), f)
            if s_name in data:
                # Two runs on the same input (say `filter` and `damage`): name this one after
                # its summary file instead, as `pmdshodh report` does.
                log.debug(f"Two summaries for input {s_name}; naming {f['fn']} as {f['s_name']}")
                s_name = f["s_name"]
            if s_name in data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, s_name)
            self.add_software_version(parsed["version"], s_name)
            data[s_name] = parsed

        data = self.ignore_samples(data)
        if len(data) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(data)} reports")

        self.write_data_file(
            {
                s: {k: v for k, v in d.items() if not isinstance(v, dict)}
                | {f"skipped_{r}": n for r, n in d["skipped"].items()}
                for s, d in data.items()
            },
            "multiqc_pmdshodh",
        )
        self.general_stats(data)
        self.read_counts_section(data)
        self.damage_section(data)
        self.score_section(data)
        self.length_section(data)

    def general_stats(self, data: dict[str, dict]):
        self.general_stats_addcols(
            {s: {k: d[k] for k in ("scored", "pct_kept", "first_c_to_t")} for s, d in data.items()},
            {
                "first_c_to_t": {
                    "title": "5′ C>T",
                    "description": "C>T rate at the first base of the 5′ end, of scored reads (pmdshodh)",
                    "suffix": "%",
                    "min": 0,
                    "max": 100,
                    "scale": "OrRd",
                    "format": "{:,.2f}",
                },
                "pct_kept": {
                    "title": "% PMD kept",
                    "description": "Scored reads kept by the PMD score threshold (pmdshodh filter)",
                    "suffix": "%",
                    "min": 0,
                    "max": 100,
                    "scale": "RdYlGn",
                },
                "scored": {
                    "title": f"PMD scored ({config.read_count_prefix})",
                    "description": f"Reads given a PMD score ({config.read_count_desc})",
                    "shared_key": "read_count",
                    "modify": lambda x: x * config.read_count_multiplier,
                    "scale": "Blues",
                    "hidden": True,
                },
            },
            namespace="pmdshodh",
        )

    def read_counts_section(self, data: dict[str, dict]):
        cats = {
            "kept": {"name": "Kept", "color": "#1a9850"},
            "removed_by_threshold": {"name": "Removed by threshold", "color": "#d73027"},
            "scored_unfiltered": {"name": "Scored (not filtered)", "color": "#4575b4"},
        }
        plot_data = {}
        for s, d in data.items():
            row = {}
            if d["kept"] is not None:
                row["kept"] = d["kept"]
                row["removed_by_threshold"] = d["removed_by_threshold"]
            else:
                row["scored_unfiltered"] = d["scored"]
            for reason, n in d["skipped"].items():
                key = f"skipped_{reason}"
                cats.setdefault(key, {"name": f"Skipped: {reason.replace('_', ' ')}"})
                row[key] = n
            plot_data[s] = row
        self.add_section(
            name="Read counts",
            anchor="pmdshodh-reads",
            description="What happened to each input read: kept or removed by the PMD score threshold, "
            "scored without filtering, or skipped (by reason).",
            helptext="Skip reasons are those recorded under `reads.skipped` in the run summary, "
            "such as unmapped reads or reads below the mapping-quality filter. Only scored reads "
            "are compared with the threshold, and only `pmdshodh filter` applies one.",
            plot=bargraph.plot(
                plot_data,
                cats,
                {
                    "id": "pmdshodh_reads_plot",
                    "title": "pmdshodh: Read counts",
                    "ylab": "Reads",
                    "cpswitch_counts_label": "Number of reads",
                },
            ),
        )

    def damage_section(self, data: dict[str, dict]):
        five = {s: d["damage"]["five_prime"] for s, d in data.items()}
        three = {s: d["damage"]["three_prime"] for s, d in data.items()}
        self.add_section(
            name="Damage by position",
            anchor="pmdshodh-damage",
            description="Substitution rate by position from each molecule end, in scored reads: "
            "5′ C>T, and 3′ G>A for double-stranded or 3′ C>T for single-stranded libraries.",
            helptext="Each rate is the count of the read base over all read bases at that reference "
            "base and position, computed from the counts in the run summary.",
            plot=linegraph.plot(
                [five, three],
                {
                    "id": "pmdshodh_damage_plot",
                    "title": "pmdshodh: Damage by position",
                    "xlab": "Position from the molecule end",
                    "ylab": "% substitution",
                    "ymin": 0,
                    "xmin": 1,
                    "tt_label": "{point.x}: {point.y:.2f}%",
                    "data_labels": [
                        {"name": "5′ C>T", "ylab": "% C>T"},
                        {"name": "3′ G>A (ds) or C>T (ss)", "ylab": "% G>A or C>T"},
                    ],
                },
            ),
        )

    def score_section(self, data: dict[str, dict]):
        hist = {s: d["score_histogram"] for s, d in data.items() if d["score_histogram"]}
        if not hist:
            return
        self.add_section(
            name="PMD score distribution",
            anchor="pmdshodh-scores",
            description="PMD scores of scored reads, in 0.5-wide bins, as a percentage of scored reads. "
            "Scores beyond the plotted range are not shown.",
            helptext="The PMD score is a log-likelihood ratio comparing a damage model with a "
            "no-damage model for each read (Skoglund et al. 2014). Reads from damaged, ancient "
            "molecules tend to score above 0.",
            plot=linegraph.plot(
                hist,
                {
                    "id": "pmdshodh_score_plot",
                    "title": "pmdshodh: PMD score distribution",
                    "xlab": "PMD score (bin centre)",
                    "ylab": "% of scored reads",
                    "ymin": 0,
                },
            ),
        )

    def length_section(self, data: dict[str, dict]):
        scored = {s: d["lengths_scored"] for s, d in data.items() if d["lengths_scored"]}
        if not scored:
            return
        kept = {s: d["lengths_kept"] for s, d in data.items() if d["lengths_kept"]}
        datasets, labels = [scored], [{"name": "Scored reads", "ylab": "Scored reads"}]
        if kept:
            datasets.append(kept)
            labels.append({"name": "Kept reads", "ylab": "Kept reads"})
        self.add_section(
            name="Read lengths",
            anchor="pmdshodh-lengths",
            description="Lengths of scored reads and, for `pmdshodh filter`, of kept reads.",
            helptext="Comparing kept with scored reads shows whether the PMD score threshold "
            "changes the length distribution.",
            plot=linegraph.plot(
                datasets,
                {
                    "id": "pmdshodh_length_plot",
                    "title": "pmdshodh: Read lengths",
                    "xlab": "Read length (bases)",
                    "ylab": "Reads",
                    "ymin": 0,
                    "data_labels": labels,
                },
            ),
        )
