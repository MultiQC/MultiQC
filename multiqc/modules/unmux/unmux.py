import logging
from collections import Counter
from typing import Any, Dict, List, Optional, Union, cast

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound, SampleGroupingConfig
from multiqc.plots import bargraph, table
from multiqc.plots.table_object import ColumnDict, InputRow
from multiqc.types import ColumnKey, LoadedFileDict, SampleGroup, SampleName

log = logging.getLogger(__name__)

Value = Union[int, float, str, None]
Metrics = Dict[str, Value]

TEXT_COLUMNS = ("pool", "sample", "sub_sample", "mean_read_len_by_input")
SUMMARY_COUNTS = ("total_reads", "reads_assigned", "reads_pass_through", "reads_unassigned", "reads_removed")
PER_SAMPLE_COLUMNS = ("pool", "sample", "sub_sample", "reads", "frac_of_pool", "frac_bases_unextracted")
SUMMARY_COLUMNS = (*SUMMARY_COUNTS, "frac_assigned", "frac_bases_unextracted", "mean_read_len_by_input")

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

    The module reads the two metrics files unmux writes when asked, recognising each by its header, so any file
    name works:

    - The `--metrics-summary` file: one row per pool, counting reads that were assigned to a sample, passed
      through, left unassigned or removed.
    - The `--metrics-per-sample` file: one row per sample, or per sub-sample, with its reads and its share of
      the pool.

    Files from many runs can be given together, and a file may hold the rows of several pools, for example when
    the files of several runs are concatenated. A pass-through run writes an empty per-sample file, which is fine.

    A compatible run:

    ```bash
    unmux r1.fq r2.fq --pool pool1 --group bc=barcodes.tsv --sample-from-group bc --out %sample.bam \\
      --metrics-per-sample %pool.unmux.per_sample.tsv \\
      --metrics-summary %pool.unmux.summary.tsv
    ```

    #### Sample names

    Pools are named by the `pool` column (unmux's `--pool`). Each sample is named by its `sample` column, as
    `<sample> (<sub_sample>)` when it has a sub-sample. If the same name turns up in more than one pool, or
    matches a pool name, the pool is added in front: `<pool>: <sample>`.

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
            info="Flexible record parsing and demultiplexing to FASTX/SAM/BAM/CRAM, splitcode-style.",
            license="MIT License",
            license_url="https://github.com/clintval/unmux/blob/main/LICENSE",
        )

        summaries = dict(sorted(self.parse_summaries().items()))
        per_sample = dict(sorted(self.parse_per_sample().items()))
        if not summaries and not per_sample:
            raise ModuleNoSamplesFound
        samples = name_samples(per_sample, set(summaries))
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

    def pool_name(self, f: LoadedFileDict[str], row: Dict[str, str]) -> str:
        return self.clean_s_name(row["pool"], f) if row.get("pool") else f["s_name"]

    def parse_summaries(self) -> Dict[str, Metrics]:
        data: Dict[str, Metrics] = {}
        for f in self.find_log_files("unmux/summary"):
            for row in read_rows(f):
                pool = self.pool_name(f, row)
                if pool in data:
                    log.debug(f"Duplicate pool found in {f['fn']}! Overwriting: {pool}")
                d = {col: to_value(col, val) for col, val in row.items() if col != "pool"}
                if any(not isinstance(d.get(col), int) for col in SUMMARY_COUNTS):
                    log.warning(f"Skipping pool '{pool}' in {f['fn']}: a read count column is missing or empty")
                    continue
                total = number(d["total_reads"])
                for col in ("reads_pass_through", "reads_unassigned", "reads_removed"):
                    d[col.replace("reads_", "frac_")] = number(d[col]) / total if total > 0 else None
                data[pool] = d
                self.add_data_source(f, pool, section="summary")
        return self.ignore_samples(data)

    def parse_per_sample(self) -> Dict[str, List[Metrics]]:
        data: Dict[str, List[Metrics]] = {}
        seen: Dict[str, str] = {}
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
                if not isinstance(d.get("reads"), int):
                    log.warning(f"Skipping sample '{row.get('sample')}' of pool '{pool}' in {f['fn']}: no reads")
                    continue
                data.setdefault(pool, []).append(d)
                self.add_data_source(f, pool, section="per_sample")
        return data

    def add_general_stats(self, summaries: Dict[str, Metrics], samples: Dict[str, Dict[str, Metrics]]) -> None:
        data: Dict[str, Dict[str, Value]] = {}
        groups: Dict[str, List[str]] = {}
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

        headers = {
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
        group_by_pool = (getattr(config, "unmux_config", None) or {}).get("general_stats_group_by_pool", True)
        self.general_stats_addcols(
            {s_name: {k: v for k, v in d.items() if v is not None} for s_name, d in data.items()},
            headers,
            group_samples_config=SampleGroupingConfig(
                explicit_groups=groups if group_by_pool else None,
                cols_to_sum=[ColumnKey("reads"), ColumnKey("frac_of_pool")],
            ),
        )

    def add_outcomes_section(self, summaries: Dict[str, Metrics]) -> None:
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

    def add_reads_per_sample_section(
        self, summaries: Dict[str, Metrics], samples: Dict[str, Dict[str, Metrics]]
    ) -> None:
        reads: Dict[str, Dict[str, float]] = {}
        percent: Dict[str, Dict[str, float]] = {}
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

    def add_metrics_section(self, summaries: Dict[str, Metrics], samples: Dict[str, Dict[str, Metrics]]) -> None:
        rows_by_group: Dict[Union[str, SampleGroup], List[InputRow]] = {}
        extra_columns: Dict[str, None] = {}
        for pool in [*summaries, *(p for p in samples if p not in summaries)]:
            by_name = samples.get(pool, {})
            rows = [InputRow(sample=SampleName(pool), data=pool_row(summaries.get(pool), list(by_name.values())))]
            for s_name, d in by_name.items():
                extras = {k: v for k, v in d.items() if k not in PER_SAMPLE_COLUMNS and v is not None}
                extra_columns.update(dict.fromkeys(extras))
                rows.append(InputRow(sample=SampleName(s_name), data={**metric_values(d), **extras}))
            rows_by_group[SampleGroup(pool)] = rows
        for d in summaries.values():
            extra_columns.update(dict.fromkeys(k for k in d if k not in SUMMARY_COLUMNS and not k.startswith("frac_")))

        any_removed = any(number(d["reads_removed"]) > 0 for d in summaries.values())
        any_pass_through = any(number(d["reads_pass_through"]) > 0 for d in summaries.values())
        headers: Dict[str, Dict[str, Any]] = {
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
                "such as barcodes, linkers and adapters",
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
            usually is the same for every sample of a pool, and a change between runs of the same assay points to
            a changed read structure.
            """,
            plot=table.plot(
                rows_by_group,
                cast(Dict[str, ColumnDict], headers),
                {
                    "id": f"{self.anchor}_metrics_table",
                    "title": "unmux: Demultiplexing Metrics",
                    "namespace": "unmux",
                },
            ),
        )


def name_samples(per_sample: Dict[str, List[Metrics]], pools: set) -> Dict[str, Dict[str, Metrics]]:
    """Name each sample `<sample> (<sub_sample>)`, adding `<pool>: ` in front when a name is not unique."""

    def base(d: Metrics) -> str:
        sample, sub_sample = str(d.get("sample") or ""), str(d.get("sub_sample") or "")
        return f"{sample} ({sub_sample})" if sub_sample else sample

    pools = pools | set(per_sample)
    counts = Counter(base(d) for rows in per_sample.values() for d in rows)
    named: Dict[str, Dict[str, Metrics]] = {}
    for pool, rows in per_sample.items():
        named[pool] = {}
        for d in rows:
            name = base(d)
            if counts[name] > 1 or name in pools:
                name = f"{pool}: {name}"
            named[pool][name] = d
    return named


def pool_row(summary: Optional[Metrics], samples: List[Metrics]) -> Dict[str, Value]:
    """The headline row of a pool: its summary, with the read columns summed over its samples."""
    if summary is None:
        fracs = [d.get("frac_of_pool") for d in samples]
        return {
            "reads": sum(number(d["reads"]) for d in samples),
            "frac_of_pool": sum(number(v) for v in fracs) if all(v is not None for v in fracs) else None,
        }
    total = number(summary["total_reads"])
    kept = number(summary["reads_assigned"]) + number(summary["reads_pass_through"])
    row: Dict[str, Value] = {
        **{k: v for k, v in summary.items() if v is not None and not k.startswith("reads_")},
        "reads": kept,
        "frac_of_pool": kept / total if total > 0 else None,
    }
    return {k: v for k, v in row.items() if v is not None}


def metric_values(d: Metrics) -> Dict[str, Value]:
    return {k: d[k] for k in ("reads", "frac_of_pool", "frac_bases_unextracted") if d.get(k) is not None}


def read_rows(f: LoadedFileDict[str]) -> List[Dict[str, str]]:
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
            log.warning(f"Skipping {f['fn']} line {line_num}: {len(fields)} fields, expected {len(header)}")
            continue
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


def number(value: Value) -> float:
    return value if isinstance(value, (int, float)) else 0
