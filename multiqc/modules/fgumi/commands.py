"""Small per-command metrics: `clip --metrics`, `filter --stats`, `copy-umi --metrics`, `retag --metrics` and
`downsample --histogram-kept/--histogram-rejected`."""

import itertools
import logging
from typing import Any, Dict, Optional, Set, Union

from pydantic import ValidationError

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph, table

from .schemas import ClippingMetric, CopyUmiMetric, DownsampleHistogramMetric, FilterStatsMetric, RetagMetric
from .util import drop_none, flatten, load_rows, pct, register, sample_name

log = logging.getLogger(__name__)

_CLIP_REASONS = {
    "five_prime": "5' end",
    "three_prime": "3' end",
    "overlapping": "Overlapping mate",
    "extending": "Extending past mate",
}


def parse_reports(module: BaseMultiqcModule) -> Set[str]:
    return _clip(module) | _filter(module) | _copy_umi(module) | _retag(module) | _downsample(module)


def _clip(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[str, Dict[str, int]]] = {}
    for f in module.find_log_files("fgumi/clip"):
        rows = load_rows(f, ClippingMetric)
        if rows is None:
            continue
        total = next((r for r in rows if r.read_type == "All"), None)
        if total is None:
            log.warning(f"Skipping {f['fn']}: fgumi clip metrics have no 'All' row")
            continue
        data[register(module, f)] = {
            f"{unit}_clipped": {reason: getattr(total, f"{unit}_clipped_{reason}") for reason in _CLIP_REASONS}
            for unit in ("reads", "bases")
        }
    data = module.ignore_samples(data)
    if not data:
        return set()
    module.add_section(
        name="Clipping",
        anchor="fgumi-clip",
        description="Reads and bases clipped by `fgumi clip` (or fgbio ClipBam), by reason (all read types combined).",
        helptext="Clipping reasons can overlap for a read, so the bars are shown side by side rather than stacked. "
        "Overlapping-mate clipping removes bases a read shares with its mate; extending clipping removes bases past "
        "the mate's end.",
        plot=bargraph.plot(
            [{s: v["reads_clipped"] for s, v in data.items()}, {s: v["bases_clipped"] for s, v in data.items()}],
            [{k: {"name": n} for k, n in _CLIP_REASONS.items()}] * 2,
            {
                "id": "fgumi_clip",
                "title": "fgumi: Clipping",
                "stacking": "group",
                "data_labels": [{"name": "Reads", "ylab": "Reads"}, {"name": "Bases", "ylab": "Bases"}],
            },
        ),
    )
    module.write_data_file(flatten(data), "multiqc_fgumi_clip")
    return set(data)


def _read_filter_stats(f: Any) -> Optional[FilterStatsMetric]:
    """The single row of a `filter --stats` file, in either layout: the headered one-row TSV, or the headerless
    key/value rows that fgumi 0.7.0 and earlier write. ``None`` (with a warning) if neither parses."""
    if f["f"] is None:
        return None
    if f["f"].split("\n", 1)[0].startswith("total_reads\tpassed_reads"):
        rows = load_rows(f, FilterStatsMetric)
        return rows[0] if rows else None
    try:
        fields = dict(line.rstrip("\r").split("\t", 1) for line in f["f"].splitlines() if line.strip())
        return FilterStatsMetric.model_validate(fields)
    except (ValueError, ValidationError) as error:
        log.warning(f"Skipping {f['fn']}: not a valid fgumi filter --stats file: {error}")
        return None


def _filter(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[str, Optional[float]]] = {}
    for f in module.find_log_files("fgumi/filter_stats"):
        row = _read_filter_stats(f)
        if row is None:
            continue
        data[register(module, f)] = {
            "passed": row.passed_reads,
            "failed": row.failed_reads,
            "pass_rate": pct(row.pass_rate),
        }
    data = module.ignore_samples(data)
    if not data:
        return set()
    module.add_section(
        name="Consensus read filtering",
        anchor="fgumi-filter",
        description="Consensus reads kept and rejected by `fgumi filter`.",
        helptext="Written by `fgumi filter --stats`. Both the headered layout and the headerless key/value layout of fgumi "
        "0.7.0 and earlier are read.",
        plot=bargraph.plot(
            {s: {"passed": v["passed"], "failed": v["failed"]} for s, v in data.items()},
            {"passed": {"name": "Passed"}, "failed": {"name": "Failed"}},
            {"id": "fgumi_filter", "title": "fgumi: Consensus read filtering", "ylab": "Reads"},
        ),
    )
    module.general_stats_addcols(
        {s: drop_none({"pass_rate": v["pass_rate"]}) for s, v in data.items()},
        {
            "pass_rate": {
                "title": "% pass filter",
                "description": "Consensus reads passing fgumi filter",
                "suffix": "%",
                "max": 100,
                "min": 0,
                "hidden": True,
            }
        },
    )
    module.write_data_file(flatten(data), "multiqc_fgumi_filter")
    return set(data)


def _copy_umi(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[str, int]] = {}
    for f in module.find_log_files("fgumi/copy_umi"):
        rows = load_rows(f, CopyUmiMetric)
        if not rows:
            continue
        row = rows[0]
        data[register(module, f)] = {
            "rx_written": row.rx_written,
            "rx_overwritten": row.rx_overwritten,
            "names_trimmed": row.names_trimmed,
        }
    data = module.ignore_samples(data)
    if not data:
        return set()
    module.add_section(
        name="Copy UMI",
        anchor="fgumi-copy-umi",
        description="Records updated by `fgumi copy-umi`.",
        helptext="Written by `fgumi copy-umi --metrics`. RX overwritten counts records that already had an RX tag; read "
        "names trimmed counts records whose UMI was removed from the name with `--remove-umi`.",
        plot=bargraph.plot(
            data,
            {
                "rx_written": {"name": "RX written"},
                "rx_overwritten": {"name": "RX overwritten"},
                "names_trimmed": {"name": "Read names trimmed"},
            },
            {"id": "fgumi_copy_umi", "title": "fgumi: Copy UMI", "ylab": "Records", "stacking": "group"},
        ),
    )
    module.write_data_file(data, "multiqc_fgumi_copy_umi")
    return set(data)


def _retag(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[str, Union[int, str]]] = {}
    samples: Set[str] = set()
    for f in module.find_log_files("fgumi/retag"):
        rows = load_rows(f, RetagMetric)
        if not rows:
            continue
        # Rows are keyed "<sample> (<operation>)", so sample filters are applied to the bare sample name here.
        s_name = sample_name(module, f)
        if module.is_ignore_sample(s_name):
            continue
        samples.add(register(module, f, s_name))
        for row in rows:
            key = f"{s_name} ({row.operation})"
            if key in data:
                # The same operation given twice: keep both rows, numbered in file order.
                key = next(k for k in (f"{s_name} ({row.operation} #{n})" for n in itertools.count(2)) if k not in data)
            data[key] = {
                "kind": row.kind,
                "records_applied": row.records_applied,
                "dst_overwritten": row.dst_overwritten,
                "src_missing": row.src_missing,
            }
    if not data:
        return set()
    module.add_section(
        name="Retag",
        anchor="fgumi-retag",
        description="Tag operations applied by `fgumi retag`, one row per sample and operation.",
        helptext="Written by `fgumi retag --metrics`: one row per `SRC::copy|move|delete::DST` operation. Source missing "
        "counts records the operation skipped because the source tag was absent.",
        plot=table.plot(
            data,
            {
                "records_applied": {"title": "Applied", "format": "{:,.0f}", "scale": "Greens"},
                "dst_overwritten": {"title": "Destination overwritten", "format": "{:,.0f}", "scale": "Oranges"},
                "src_missing": {"title": "Source missing", "format": "{:,.0f}", "scale": "Reds"},
                "kind": {"title": "Operation"},
            },
            {"id": "fgumi_retag", "title": "fgumi: Retag", "namespace": "fgumi retag"},
        ),
    )
    module.write_data_file(data, "multiqc_fgumi_retag")
    return samples


def _downsample(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[int, int]] = {}
    for f in module.find_log_files("fgumi/downsample_histogram"):
        rows = load_rows(f, DownsampleHistogramMetric)
        if rows is None:
            continue
        data[register(module, f)] = {r.family_size: r.count for r in rows}
    data = module.ignore_samples(data)
    if not data:
        return set()
    module.add_section(
        name="Downsampled family sizes",
        anchor="fgumi-downsample",
        description="Family sizes kept or rejected by `fgumi downsample`, one line per histogram file.",
        helptext="Written by `fgumi downsample --histogram-kept` / `--histogram-rejected`. Each file is its own line, "
        "named after the file.",
        plot=linegraph.plot(
            data,
            {
                "id": "fgumi_downsample",
                "title": "fgumi: Downsampled family sizes",
                "xlab": "Family size",
                "ylab": "Families",
                "xlog": True,
            },
        ),
    )
    module.write_data_file(data, "multiqc_fgumi_downsample")
    return set(data)
