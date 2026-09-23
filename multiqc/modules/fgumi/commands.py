"""Small per-command metrics: `clip --metrics`, `filter --stats`, `copy-umi --metrics`, `retag --metrics` and
`downsample --histogram-kept/--histogram-rejected`."""

from typing import Dict, Set

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph, table

from .schemas import ClippingMetric, CopyUmiMetric, DownsampleHistogramMetric, FilterStatsMetric, RetagMetric
from .util import load_rows, pct, sample_name

_CLIP_REASONS = {
    "five_prime": "5' end",
    "three_prime": "3' end",
    "overlapping": "Overlapping mate",
    "extending": "Extending past mate",
}


def parse_reports(module: BaseMultiqcModule) -> Set[str]:
    return _clip(module) | _filter(module) | _copy_umi(module) | _retag(module) | _downsample(module)


def _register(module: BaseMultiqcModule, f: Dict) -> str:
    s_name = sample_name(module, f)
    module.add_data_source(f, s_name)
    module.add_software_version(None, s_name)
    return s_name


def _clip(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[str, Dict[str, int]]] = {}
    for f in module.find_log_files("fgumi/clip"):
        rows = load_rows(f, ClippingMetric)
        total = next((r for r in rows or [] if r.read_type == "All"), None)
        if total is None:
            continue
        data[_register(module, f)] = {
            unit: {reason: getattr(total, f"{unit}_clipped_{reason}") for reason in _CLIP_REASONS}
            for unit in ("reads", "bases")
        }
    data = module.ignore_samples(data)
    if not data:
        return set()
    module.add_section(
        name="Clipping",
        anchor="fgumi-clip",
        description="Reads and bases clipped by `fgumi clip`, by reason (all read types combined).",
        plot=bargraph.plot(
            [{s: v["reads"] for s, v in data.items()}, {s: v["bases"] for s, v in data.items()}],
            [{k: {"name": n} for k, n in _CLIP_REASONS.items()}] * 2,
            {
                "id": "fgumi_clip",
                "title": "fgumi: Clipping",
                "stacking": "group",
                "data_labels": [{"name": "Reads", "ylab": "Reads"}, {"name": "Bases", "ylab": "Bases"}],
            },
        ),
    )
    module.write_data_file(data, "multiqc_fgumi_clip")
    return set(data)


def _filter(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[str, object]] = {}
    for f in module.find_log_files("fgumi/filter_stats"):
        rows = load_rows(f, FilterStatsMetric)
        if not rows:
            continue
        row = rows[0]
        data[_register(module, f)] = {
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
        plot=bargraph.plot(
            {s: {"passed": v["passed"], "failed": v["failed"]} for s, v in data.items()},
            {"passed": {"name": "Passed"}, "failed": {"name": "Failed"}},
            {"id": "fgumi_filter", "title": "fgumi: Consensus read filtering", "ylab": "Reads"},
        ),
    )
    module.general_stats_addcols(
        {s: {"pass_rate": v["pass_rate"]} for s, v in data.items()},
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
    module.write_data_file(data, "multiqc_fgumi_filter")
    return set(data)


def _copy_umi(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[str, int]] = {}
    for f in module.find_log_files("fgumi/copy_umi"):
        rows = load_rows(f, CopyUmiMetric)
        if not rows:
            continue
        row = rows[0]
        data[_register(module, f)] = {
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
    data: Dict[str, Dict[str, object]] = {}
    for f in module.find_log_files("fgumi/retag"):
        rows = load_rows(f, RetagMetric)
        if not rows:
            continue
        s_name = _register(module, f)
        for row in rows:
            data[f"{s_name} ({row.operation})"] = {
                "kind": row.kind,
                "records_applied": row.records_applied,
                "dst_overwritten": row.dst_overwritten,
                "src_missing": row.src_missing,
            }
    data = module.ignore_samples(data)
    if not data:
        return set()
    module.add_section(
        name="Retag",
        anchor="fgumi-retag",
        description="Tag operations applied by `fgumi retag`, one row per sample and operation.",
        plot=table.plot(
            data,
            {
                "records_applied": {"title": "Applied", "format": "{:,.0f}"},
                "dst_overwritten": {"title": "Destination overwritten", "format": "{:,.0f}"},
                "src_missing": {"title": "Source missing", "format": "{:,.0f}"},
                "kind": {"title": "Operation"},
            },
            {"id": "fgumi_retag", "title": "fgumi: Retag", "namespace": "fgumi retag"},
        ),
    )
    module.write_data_file(data, "multiqc_fgumi_retag")
    return {key.rsplit(" (", 1)[0] for key in data}


def _downsample(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[int, int]] = {}
    for f in module.find_log_files("fgumi/downsample_histogram"):
        rows = load_rows(f, DownsampleHistogramMetric)
        if rows is None:
            continue
        data[_register(module, f)] = {r.family_size: r.count for r in rows}
    data = module.ignore_samples(data)
    if not data:
        return set()
    module.add_section(
        name="Downsampled family sizes",
        anchor="fgumi-downsample",
        description="Family sizes kept or rejected by `fgumi downsample`, one line per histogram file.",
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
