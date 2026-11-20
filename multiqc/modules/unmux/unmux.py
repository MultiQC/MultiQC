import logging
from collections import Counter
from collections.abc import Callable, Sequence
from typing import Any, cast

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound, SampleGroupingConfig
from multiqc.plots import bargraph, table
from multiqc.plots.table_object import ColumnDict, InputRow
from multiqc.types import ColumnKey, LoadedFileDict, SampleGroup, SampleName

log = logging.getLogger(__name__)

Value = int | float | str | None
Metrics = dict[str, Value]

TEXT_COLUMNS = ("pool", "sample", "sub_sample", "mean_read_len_by_input", "reason", "group", "stream")
SUMMARY_COUNTS = ("total_reads", "reads_assigned", "reads_pass_through", "reads_unassigned", "reads_removed")
PER_SAMPLE_COLUMNS = ("pool", "sample", "sub_sample", "reads", "frac_of_pool", "frac_bases_unextracted")
SUMMARY_COLUMNS = (*SUMMARY_COUNTS, "frac_assigned", "frac_bases_unextracted", "mean_read_len_by_input")
# Half of the 0.1% the table shows, so a sample's unextracted share is shown only where it visibly differs.
UNEXTRACTED_TOLERANCE = 0.0005

REASONS = {
    "min_finds_per_group": "too few matches",
    "max_finds_per_group": "too many matches",
    "min_finds_per_tag": "a tag found too few times",
    "max_finds_per_tag": "a tag found too many times",
}
NO_SAMPLE_COLOR = "#999999"

OUTCOME_CATS = {
    "reads_assigned": {"name": "Assigned", "color": "#437bb1"},
    "reads_pass_through": {"name": "Passed through", "color": "#7cb5ec"},
    "reads_removed": {"name": "Removed", "color": "#f7a35c"},
    "reads_unassigned": {"name": "Unassigned", "color": "#999999"},
}


class MultiqcModule(BaseMultiqcModule):
    """
    unmux parses and demultiplexes FASTX, SAM, BAM and CRAM records with splitcode-style tag groups. One run
    reads one pool and splits it into samples, and optionally each sample into sub-samples (written as the read
    group `LB`). A run with no samples passes every record through, and `--remove` rules drop records that match
    a tag.

    The module reads the three metrics files unmux writes when asked, recognising each by its header, so any file
    name works:

    - The `--metrics-summary` file: one row per pool, counting reads that were assigned to a sample, passed
      through, left unassigned or removed.
    - The `--metrics-per-sample` file: one row per sample, or per sub-sample, with its reads and its share of
      the pool.
    - The `--metrics-unassigned` file: one row per reason the pool's unassigned reads matched no sample, such
      as a tag group with too few matches.

    Files from many runs can be given together, and a file may hold the rows of several pools, for example when
    the files of several runs are concatenated. A pass-through run writes an empty per-sample file, which is fine.

    A compatible run:

    ```bash
    unmux r1.fq r2.fq --pool pool1 --group bc=barcodes.tsv --sample-from-group bc --out %sample.bam \\
      --metrics-per-sample %pool.unmux.per_sample.tsv \\
      --metrics-summary %pool.unmux.summary.tsv \\
      --metrics-unassigned %pool.unmux.unassigned.tsv
    ```

    #### Sample names

    Pools are named by the `pool` column (unmux's `--pool`). Each sample is named by its `sample` column, as
    `<sample> (<sub_sample>)` when it has a sub-sample. If the same name turns up in more than one pool, or
    matches a pool name, the pool is added in front: `<pool>: <sample>`.

    Sample names are cleaned like any other MultiQC sample name, so `sample_names_replace` applies to them. For
    example, to drop the sub-sample so that unmux's rows join the rows other modules write for each sample:

    ```yaml
    sample_names_replace_regex: true
    sample_names_replace:
      ' \\([^)]*\\)$': ""
    unmux_config:
      general_stats_group_by_pool: false
    ```

    #### Grouping by pool

    The tables and the reads per sample plot nest each sample under its pool. In the General Statistics table,
    a pool's samples are shown as one collapsed row named after the pool, which you can click to expand. As
    MultiQC joins General Statistics rows by name, this grouping keeps unmux's sample rows apart from rows that
    other modules write for the same samples. To show one flat row per sample instead, use:

    ```yaml
    unmux_config:
      general_stats_group_by_pool: false
    ```
    """

    def __init__(self):
        super().__init__(
            name="unmux",
            anchor="unmux",
            href="https://github.com/clintval/unmux",
            info="Flexible record parsing and demultiplexing to FASTX/SAM/BAM/CRAM.",
            # doi="", No DOI available
            license="MIT License",
            license_url="https://github.com/clintval/unmux/blob/main/LICENSE",
        )

        summaries = dict(sorted(self.parse_summaries().items()))
        self.files_by_row: dict[int, LoadedFileDict[str]] = {}
        per_sample = dict(sorted(self.parse_per_sample().items()))
        unassigned = dict(sorted(self.parse_unassigned().items()))
        if not summaries and not per_sample and not unassigned:
            raise ModuleNoSamplesFound
        samples = name_samples(
            per_sample, set(summaries), clean=lambda name, d: self.clean_s_name(name, self.files_by_row[id(d)])
        )
        samples = {
            pool: {s_name: d for s_name, d in by_name.items() if not self.is_ignore_sample(s_name)}
            for pool, by_name in samples.items()
        }
        log.info(
            f"Found {len(summaries)} pool summaries and {sum(len(s) for s in samples.values())} samples "
            f"in {len(samples)} pools"
        )

        self.add_general_stats(summaries, samples)
        if summaries:
            self.add_outcomes_section(summaries)
        if any(unassigned.values()):
            self.add_unassigned_section(unassigned)
        if any(samples.values()):
            self.add_reads_per_sample_section(summaries, samples)
        self.add_metrics_section(summaries, samples)

        if summaries:
            self.write_data_file(summaries, "multiqc_unmux_summary")
        if samples:
            self.write_data_file(
                {s_name: d for by_name in samples.values() for s_name, d in by_name.items()},
                "multiqc_unmux_per_sample",
            )
        if unassigned:
            self.write_data_file(
                {pool: {reason_label(d): d["reads"] for d in rows} for pool, rows in unassigned.items()},
                "multiqc_unmux_unassigned",
            )

    def pool_name(self, f: LoadedFileDict[str], row: dict[str, str]) -> str:
        return self.clean_s_name(row["pool"], f) if row.get("pool") else f["s_name"]

    def parse_summaries(self) -> dict[str, Metrics]:
        data: dict[str, Metrics] = {}
        for f in self.find_log_files("unmux/summary"):
            for row in read_rows(f):
                pool = self.pool_name(f, row)
                if pool in data:
                    log.debug(f"Duplicate pool found in {f['fn']}! Overwriting: {pool}")
                d = {col: to_value(col, val) for col, val in row.items() if col != "pool"}
                require_counts(f, d, SUMMARY_COUNTS)
                total = number(d["total_reads"])
                for col in ("reads_pass_through", "reads_unassigned", "reads_removed"):
                    d[col.replace("reads_", "frac_")] = number(d[col]) / total if total > 0 else None
                data[pool] = d
                self.add_data_source(f, pool, section="summary")
                # unmux does not write its version into the metrics files.
                self.add_software_version(None, pool)
        return self.ignore_samples(data)

    def parse_per_sample(self) -> dict[str, list[Metrics]]:
        data: dict[str, list[Metrics]] = {}
        seen: dict[str, str] = {}
        for f in self.find_log_files("unmux/per_sample"):
            for row in read_rows(f):
                pool = self.pool_name(f, row)
                if self.is_ignore_sample(pool):
                    continue
                # A later file for the same pool replaces the earlier one rather than adding to it.
                if seen.setdefault(pool, f["fn"]) != f["fn"]:
                    log.debug(f"Duplicate pool found in {f['fn']}! Overwriting: {pool}")
                    data[pool] = []
                    seen[pool] = f["fn"]
                d = {col: to_value(col, val) for col, val in row.items()}
                d["pool"] = pool
                require_counts(f, d, ("reads",))
                self.files_by_row[id(d)] = f
                data.setdefault(pool, []).append(d)
                self.add_data_source(f, pool, section="per_sample")
                self.add_software_version(None, pool)
        return data

    def parse_unassigned(self) -> dict[str, list[Metrics]]:
        data: dict[str, list[Metrics]] = {}
        seen: dict[str, str] = {}
        for f in self.find_log_files("unmux/unassigned"):
            for row in read_rows(f):
                pool = self.pool_name(f, row)
                if self.is_ignore_sample(pool):
                    continue
                if seen.setdefault(pool, f["fn"]) != f["fn"]:
                    log.debug(f"Duplicate pool found in {f['fn']}! Overwriting: {pool}")
                    data[pool] = []
                    seen[pool] = f["fn"]
                d = {col: to_value(col, val) for col, val in row.items() if col != "pool"}
                require_counts(f, d, ("reads",))
                data.setdefault(pool, []).append(d)
                self.add_data_source(f, pool, section="unassigned")
                self.add_software_version(None, pool)
        return data

    def add_general_stats(self, summaries: dict[str, Metrics], samples: dict[str, dict[str, Metrics]]) -> None:
        data: dict[str, dict[str, Value]] = {}
        groups: dict[str, list[str]] = {}
        for pool, by_name in samples.items():
            for s_name, d in by_name.items():
                data[s_name] = {"reads": d["reads"], "frac_of_pool": d.get("frac_of_pool")}
            groups[pool] = list(by_name)
        for pool, d in summaries.items():
            if samples.get(pool):
                continue
            total = number(d["total_reads"])
            kept = number(d["reads_assigned"]) + number(d["reads_pass_through"])
            data[pool] = {"reads": kept, "frac_of_pool": kept / total if total > 0 else None}
        if not data:
            return

        headers: dict[str, dict[str, Any]] = {
            "frac_of_pool": {
                "title": "% of Pool",
                "description": "Share of the pool's reads written to this sample. A collapsed pool row adds up "
                "its samples, so it shows the share of the pool assigned to any sample",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "YlGn",
                "modify": lambda x: x * 100.0,
            },
            "reads": {
                "title": f"{config.read_count_prefix} Reads",
                "description": f"Reads written to this sample ({config.read_count_desc})",
                "min": 0,
                "format": "{:,.2f}",
                "scale": "Blues",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
                "hidden": True,
            },
        }
        general_stats_headers = self.get_general_stats_headers(all_headers=cast(dict[str, ColumnDict], headers))
        if not general_stats_headers:
            return
        group_by_pool = (getattr(config, "unmux_config", None) or {}).get("general_stats_group_by_pool", True)
        self.general_stats_addcols(
            {s_name: {k: v for k, v in d.items() if v is not None} for s_name, d in data.items()},
            general_stats_headers,
            group_samples_config=SampleGroupingConfig(
                explicit_groups=groups if group_by_pool else None,
                cols_to_sum=[ColumnKey("reads"), ColumnKey("frac_of_pool")],
            ),
        )

    def add_outcomes_section(self, summaries: dict[str, Metrics]) -> None:
        plot_data = {pool: {cat: number(d[cat]) for cat in OUTCOME_CATS} for pool, d in summaries.items()}
        self.add_section(
            name="Read Outcomes",
            anchor=f"{self.anchor}-outcomes",
            description="What happened to each read of each pool.",
            helptext="""
            Each read of a pool has exactly one outcome, so a pool's bars add up to its reads.

            * _Assigned_: matched a sample, or a sub-sample, and was written to its output.
            * _Passed through_: written out unchanged by a run with no samples.
            * _Removed_: matched a `--remove` rule and was dropped, or written to that rule's own file.
            * _Unassigned_: matched no sample. These go to the `--unassigned` output when one is given. A high
              share points to a wrong barcode list, a wrong tag location, or too few allowed mismatches.
            """,
            plot=bargraph.plot(
                plot_data,
                OUTCOME_CATS,
                {
                    "id": f"{self.anchor}_outcomes",
                    "title": "unmux: Read Outcomes",
                    "ylab": "Reads",
                    "cpswitch_counts_label": "Number of Reads",
                    "cpswitch_c_active": False,
                },
            ),
        )

    def add_unassigned_section(self, unassigned: dict[str, list[Metrics]]) -> None:
        reads: dict[str, dict[str, float]] = {}
        of_unassigned: dict[str, dict[str, float]] = {}
        of_pool: dict[str, dict[str, float]] = {}
        cats: dict[str, dict[str, str]] = {}
        for pool, rows in unassigned.items():
            total = sum(number(d["reads"]) for d in rows)
            for d in rows:
                label = reason_label(d)
                cats.setdefault(label, {"name": label})
                reads.setdefault(pool, {})[label] = number(d["reads"])
                if total > 0:
                    of_unassigned.setdefault(pool, {})[label] = number(d["reads"]) / total * 100.0
                of_pool.setdefault(pool, {})[label] = number(d["frac_of_pool"]) * 100.0
        # Groups first, in the order unmux matched them, then missing streams, then reads no sample claimed.
        cats = dict(sorted(cats.items(), key=lambda item: label_rank(item[0])))
        if "No sample" in cats:
            cats["No sample"]["color"] = NO_SAMPLE_COLOR

        self.add_section(
            name="Unassigned Reasons",
            anchor=f"{self.anchor}-unassigned",
            description="Why the unassigned reads of each pool matched no sample, by tag group.",
            helptext="""
            unmux matches the tag groups of a read in order and stops at the first that fails, so each unassigned
            read is counted once, under the group that failed it. A group fails when it has:

            * _too few matches_: fewer of its tags matched than `minFindsPerGroup` allows. With
              `minFindsPerGroup=1` no tag of the group matched within the allowed distance and location.
            * _too many matches_: more of its tags matched than `maxFindsPerGroup` allows, for example two
              different barcodes in one window.
            * _a tag found too few times_ or _too many times_: one tag matched fewer or more times than
              `minFindsPerTag` or `maxFindsPerTag` allow.

            A read can also pass every group and still be unassigned:

            * _stream missing_: a `--template` stream is built from a group that did not match.
            * _No sample_: every group matched, but no sample claims the tags it matched, for example a barcode
              left out of the sample sheet.

            _% of Unassigned_ divides each by the pool's unassigned reads, so a pool's bar adds up to 100%.
            _% of Pool_ divides each by all reads of the pool. A group that fails more reads in one pool than in
            the others points to a problem with that part of the read, such as a worn primer or a shifted
            position.
            """,
            plot=bargraph.plot(
                [of_unassigned, reads, of_pool],
                [cats, cats, cats],
                {
                    "id": f"{self.anchor}_unassigned",
                    "title": "unmux: Unassigned Reasons",
                    "ylab": "% of Unassigned",
                    "cpswitch": False,
                    "data_labels": [
                        {"name": "% of Unassigned", "ylab": "% of Unassigned", "ysuffix": "%", "ymax": 100},
                        {"name": "Reads", "ylab": "Reads"},
                        {"name": "% of Pool", "ylab": "% of Pool", "ysuffix": "%"},
                    ],
                },
            ),
        )

    def add_reads_per_sample_section(
        self, summaries: dict[str, Metrics], samples: dict[str, dict[str, Metrics]]
    ) -> None:
        reads: dict[str, dict[str, float]] = {}
        percent: dict[str, dict[str, float]] = {}
        for pool, by_name in samples.items():
            if not by_name:
                continue
            # The pool's unassigned and removed bars come first, so they head the pool's block of bars.
            summary = summaries.get(pool)
            total = number(summary["total_reads"]) if summary else 0
            for cat in ("reads_unassigned", "reads_removed"):
                if summary is None or number(summary[cat]) == 0:
                    continue
                key = f"{pool}: {OUTCOME_CATS[cat]['name'].lower()}"
                reads[key] = {cat: number(summary[cat])}
                if total > 0:
                    percent[key] = {cat: number(summary[cat]) / total * 100.0}
            for s_name, d in by_name.items():
                reads[s_name] = {"reads": number(d["reads"])}
                frac = d.get("frac_of_pool")
                if isinstance(frac, (int, float)):
                    percent[s_name] = {"reads": frac * 100.0}

        cats = {
            "reads": OUTCOME_CATS["reads_assigned"],
            "reads_removed": OUTCOME_CATS["reads_removed"],
            "reads_unassigned": OUTCOME_CATS["reads_unassigned"],
        }
        self.add_section(
            name="Reads per Sample",
            anchor=f"{self.anchor}-reads-per-sample",
            description="Reads written to each sample, grouped by pool, with the pool's unassigned and removed "
            "reads alongside.",
            helptext="""
            Each bar is a sample, or a sub-sample. The bars of one pool are drawn together, headed by the pool's
            unassigned and removed reads. _% of Pool_ divides each by all reads of the pool, so a pool's bars add
            up to 100%. A sample far below its neighbours
            may have failed, been under-loaded, or have a barcode that does not match the one given to unmux.
            """,
            plot=bargraph.plot(
                [reads, percent],
                [cats, cats],
                {
                    "id": f"{self.anchor}_reads_per_sample",
                    "title": "unmux: Reads per Sample",
                    "ylab": "Reads",
                    "cpswitch": False,
                    "sort_samples": False,
                    "data_labels": [
                        {"name": "Reads", "ylab": "Reads"},
                        {"name": "% of Pool", "ylab": "% of Pool", "ysuffix": "%", "ymax": 100},
                    ],
                },
            ),
        )

    def add_metrics_section(self, summaries: dict[str, Metrics], samples: dict[str, dict[str, Metrics]]) -> None:
        rows_by_group: dict[str | SampleGroup, list[InputRow]] = {}
        extra_columns: dict[str, None] = {}
        for pool in [*summaries, *(p for p in samples if p not in summaries)]:
            by_name = samples.get(pool, {})
            pool_unextracted = (summaries.get(pool) or {}).get("frac_bases_unextracted")
            rows = [InputRow(sample=SampleName(pool), data=pool_row(summaries.get(pool), list(by_name.values())))]
            for s_name, d in by_name.items():
                extras = {k: v for k, v in d.items() if k not in PER_SAMPLE_COLUMNS and v is not None}
                extra_columns.update(dict.fromkeys(extras))
                rows.append(InputRow(sample=SampleName(s_name), data={**metric_values(d, pool_unextracted), **extras}))
            rows_by_group[SampleGroup(pool)] = rows
        for d in summaries.values():
            extra_columns.update(dict.fromkeys(k for k in d if k not in SUMMARY_COLUMNS and not k.startswith("frac_")))

        any_removed = any(number(d["reads_removed"]) > 0 for d in summaries.values())
        any_pass_through = any(number(d["reads_pass_through"]) > 0 for d in summaries.values())
        headers: dict[str, dict[str, Any]] = {
            "total_reads": {
                "title": f"{config.read_count_prefix} Pool Reads",
                "description": f"All reads of the pool, shown on pool rows ({config.read_count_desc})",
                "min": 0,
                "format": "{:,.2f}",
                "scale": "Purples",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "reads": {
                "title": f"{config.read_count_prefix} Reads",
                "description": "Reads written to the sample. On a pool row, reads written to any sample or passed "
                f"through ({config.read_count_desc})",
                "min": 0,
                "format": "{:,.2f}",
                "scale": "Blues",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
            },
            "frac_of_pool": {
                "title": "% of Pool",
                "description": "Reads written to the sample as a share of all reads of the pool. On a pool row, the "
                "share written to any sample or passed through",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "YlGn",
                "modify": lambda x: x * 100.0,
            },
            "frac_unassigned": {
                "title": "% Unassigned",
                "description": "Reads that matched no sample, as a share of the pool",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
                "modify": lambda x: x * 100.0,
            },
            "frac_removed": {
                "title": "% Removed",
                "description": "Reads dropped by a --remove rule, as a share of the pool",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Oranges",
                "modify": lambda x: x * 100.0,
                "hidden": not any_removed,
            },
            "frac_pass_through": {
                "title": "% Passed Through",
                "description": "Reads written out unchanged by a run with no samples, as a share of the pool",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "PuBu",
                "modify": lambda x: x * 100.0,
                "hidden": not any_pass_through,
            },
            "frac_bases_unextracted": {
                "title": "% Unextracted",
                "description": "Bases of the kept reads that were not written into a --template or --tag output, "
                "such as barcodes, linkers and adapters. Shown on a sample row only where it differs from the pool's",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Greys",
                "modify": lambda x: x * 100.0,
            },
            "mean_read_len_by_input": {
                "title": "Read Lengths",
                "description": "Mean length of the reads of each input file, in input order",
                "scale": False,
                "hidden": True,
            },
        }
        for col in extra_columns:
            headers[col] = {"title": col, "description": f"{col}, as written by unmux", "hidden": True}

        self.add_section(
            name="Demultiplexing Metrics",
            anchor=f"{self.anchor}-metrics",
            description="Each pool, with its samples nested under it.",
            helptext="""
            Click a pool row to show its samples. Sample rows share the read columns of the pool row: a pool row's
            _Reads_ and _% of Pool_ are the sums over its samples, or over its passed through reads when the run
            had no samples. _Pool Reads_, _% Unassigned_ and _% Removed_ belong to the pool only.

            _% Unextracted_ is the share of bases left out of every `--template` and `--tag` output, such as
            barcodes, linkers and adapters. It is 0 when reads are written out whole. Set by the read layout, it
            usually is the same for every sample of a pool, so a sample row shows it only where it differs from
            its pool's, for example when that sample's reads are shorter. A change between runs of the same assay
            points to a changed read structure.
            """,
            plot=table.plot(
                rows_by_group,
                cast(dict[str, ColumnDict], headers),
                {
                    "id": f"{self.anchor}_metrics_table",
                    "title": "unmux: Demultiplexing Metrics",
                    "namespace": "unmux",
                },
            ),
        )


def reason_label(d: Metrics) -> str:
    """A category name for an unassigned reason, such as `round1: too few matches` or `No sample`."""
    reason = str(d["reason"])
    if reason == "no_sample":
        return "No sample"
    if reason == "missing_stream":
        return f"{d.get('stream') or '?'}: stream missing"
    label = REASONS.get(reason, reason)
    return f"{d['group']}: {label}" if d.get("group") else label


def label_rank(label: str) -> int:
    """Sort missing streams after the groups and reads no sample claimed last, keeping the order seen otherwise."""
    if label == "No sample":
        return 2
    return 1 if label.endswith(": stream missing") else 0


def name_samples(
    per_sample: dict[str, list[Metrics]],
    pools: set,
    clean: Callable[[str, Metrics], str] = lambda name, _: name,
) -> dict[str, dict[str, Metrics]]:
    """Name each sample `<sample> (<sub_sample>)`, cleaned, adding `<pool>: ` in front when a name is not unique."""

    def base(d: Metrics) -> str:
        sample, sub_sample = str(d.get("sample") or ""), str(d.get("sub_sample") or "")
        return clean(f"{sample} ({sub_sample})" if sub_sample else sample, d)

    pools = pools | set(per_sample)
    bases = {id(d): base(d) for rows in per_sample.values() for d in rows}
    counts = Counter(bases.values())
    named: dict[str, dict[str, Metrics]] = {}
    for pool, rows in per_sample.items():
        named[pool] = {}
        for d in rows:
            name = bases[id(d)]
            if counts[name] > 1 or name in pools:
                name = f"{pool}: {name}"
            if name in named[pool]:
                log.warning(f"Two samples of pool {pool} are both named {name!r}; keeping the last")
            named[pool][name] = d
    return named


def pool_row(summary: Metrics | None, samples: list[Metrics]) -> dict[str, Value]:
    """The headline row of a pool: its summary, with the read columns summed over its samples."""
    if summary is None:
        fracs = [d.get("frac_of_pool") for d in samples]
        return {
            "reads": sum(number(d["reads"]) for d in samples),
            "frac_of_pool": sum(number(v) for v in fracs) if all(v is not None for v in fracs) else None,
        }
    total = number(summary["total_reads"])
    kept = number(summary["reads_assigned"]) + number(summary["reads_pass_through"])
    row: dict[str, Value] = {
        **{k: v for k, v in summary.items() if v is not None and not k.startswith("reads_")},
        "reads": kept,
        "frac_of_pool": kept / total if total > 0 else None,
    }
    return {k: v for k, v in row.items() if v is not None}


def metric_values(d: Metrics, pool_unextracted: Value = None) -> dict[str, Value]:
    """A sample's read columns, with its unextracted share only where it differs from its pool's."""
    values = {k: d[k] for k in ("reads", "frac_of_pool") if d.get(k) is not None}
    unextracted = d.get("frac_bases_unextracted")
    if unextracted is None:
        return values
    if (
        isinstance(pool_unextracted, (int, float))
        and abs(number(unextracted) - pool_unextracted) < UNEXTRACTED_TOLERANCE
    ):
        return values
    return {**values, "frac_bases_unextracted": unextracted}


def read_rows(f: LoadedFileDict[str]) -> list[dict[str, str]]:
    """The rows of an unmux metrics TSV as `{column: value}`, skipping repeated headers from concatenated files."""
    lines = f["f"].splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    rows = []
    for line_num, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if not line.strip() or fields == header:
            continue
        if len(fields) != len(header):
            raise ValueError(f"unmux: {f['fn']} line {line_num} has {len(fields)} fields, expected {len(header)}")
        rows.append(dict(zip(header, fields)))
    return rows


def to_value(column: str, value: str) -> Value:
    """A cell as a number, or as text for the text columns and any column that is not numeric."""
    if value == "":
        return None
    if column in TEXT_COLUMNS:
        return value
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def require_counts(f: LoadedFileDict[str], d: Metrics, columns: Sequence[str]) -> None:
    """Fail on a read count unmux always writes but that is missing or not a whole number."""
    for column in columns:
        if not isinstance(d.get(column), int):
            raise ValueError(  # noqa: TRY004  a malformed file is a bad value, not a caller's type error
                f"unmux: {f['fn']} has no whole number in the required column '{column}'"
            )


def number(value: Value) -> float:
    if isinstance(value, (int, float)):
        return value
    raise ValueError(f"unmux: expected a number, found {value!r}")
