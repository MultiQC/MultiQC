"""Consensus family sizes from `simplex-metrics` / `duplex-metrics` (and the consensus callers' `--metrics`):
per-strand family sizes (`<prefix>.family_sizes.txt`) and the duplex AB x BA family sizes
(`<prefix>.duplex_family_sizes.txt`), drawn like fgbio CollectDuplexSeqMetrics' plots #1-#3."""

import logging
import math
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import heatmap, linegraph

from .schemas import DuplexFamilySizeMetric, DuplexStrandFamilySizeMetric, SimplexFamilySizeMetric
from .util import drop_none, load_rows, pct, safe_id, sample_name

log = logging.getLogger(__name__)

HEATMAP_CAP = 20  # sizes at or above this fold into one "20+" row/column


def parse_reports(module: BaseMultiqcModule) -> Set[str]:
    return _strand_family_sizes(module) | _ab_ba_family_sizes(module)


def _strand_family_sizes(module: BaseMultiqcModule) -> Set[str]:
    by_kind: Dict[str, Dict[str, Dict[str, Dict[int, Optional[float]]]]] = {"simplex": {}, "duplex": {}}
    for f in module.find_log_files("fgumi/strand_family_sizes"):
        if f["f"] is None:
            continue
        is_duplex = "ds_count" in f["f"].split("\n", 1)[0].split("\t")
        rows = load_rows(f, DuplexStrandFamilySizeMetric if is_duplex else SimplexFamilySizeMetric)
        if rows is None:
            continue
        s_name = sample_name(module, f)
        module.add_data_source(f, s_name)
        module.add_software_version(None, s_name)
        strands = ["cs", "ss", "ds"] if is_duplex else ["cs", "ss"]
        entry: Dict[str, Dict[int, Optional[float]]] = {}
        for strand in strands:
            entry[strand] = {r.family_size: getattr(r, f"{strand}_count") for r in rows}
            entry[f"{strand}_cumulative"] = {
                r.family_size: pct(getattr(r, f"{strand}_fraction_gt_or_eq_size")) for r in rows
            }
        by_kind["duplex" if is_duplex else "simplex"][s_name] = entry

    parsed: Set[str] = set()
    labels = {"cs": "By coordinate + strand", "ss": "Single-strand families", "ds": "Double-strand families"}
    for kind, data in by_kind.items():
        data = module.ignore_samples(data)
        if not data:
            continue
        strands = ["cs", "ss", "ds"] if kind == "duplex" else ["cs", "ss"]
        counts = [{s: drop_none(data[s][strand]) for s in data} for strand in strands]
        cumulative = [{s: drop_none(data[s][f"{strand}_cumulative"]) for s in data} for strand in strands]
        module.add_section(
            name=f"{kind.capitalize()} family sizes",
            anchor=f"fgumi-{kind}-family-sizes",
            description=f"Family size distributions from `fgumi {kind}-metrics`, one tab per grouping.",
            helptext="CS: grouped by coordinate and strand only. SS: also by UMI (single-strand families). "
            "DS: single-strand families paired into double-strand families. Matches fgbio CollectDuplexSeqMetrics.",
            plot=linegraph.plot(
                counts + cumulative,
                {
                    "id": f"fgumi_{kind}_family_sizes",
                    "title": f"fgumi: {kind.capitalize()} family sizes",
                    "xlab": "Family size (read pairs)",
                    "xlog": True,
                    "data_labels": [{"name": labels[s], "ylab": "Families"} for s in strands]
                    + [
                        {"name": f"{labels[s]} (cumulative)", "ylab": "% at or above size", "ysuffix": "%"}
                        for s in strands
                    ],
                },
            ),
        )
        module.write_data_file(
            {s: {strand: data[s][strand] for strand in strands} for s in data},
            f"multiqc_fgumi_{kind}_family_sizes",
        )
        parsed |= set(data)
    return parsed


def _label(size: int) -> str:
    return f"{HEATMAP_CAP}+" if size >= HEATMAP_CAP else str(size)


def _heatmap_matrix(
    rows: List[DuplexFamilySizeMetric], min_ba: int
) -> Optional[Tuple[List[List[Optional[float]]], List[str], List[str]]]:
    """log10(count) by (capped AB size, capped BA size); rows with ``ba_size < min_ba`` are dropped."""
    cells: Dict[Tuple[int, int], int] = defaultdict(int)
    for r in rows:
        if r.ba_size >= min_ba:
            cells[(min(r.ab_size, HEATMAP_CAP), min(r.ba_size, HEATMAP_CAP))] += r.count
    if not cells:
        return None
    ab_sizes = sorted({ab for ab, _ in cells})
    ba_sizes = sorted({ba for _, ba in cells})
    matrix = [[math.log10(cells[(ab, ba)]) if cells.get((ab, ba)) else None for ba in ba_sizes] for ab in ab_sizes]
    return matrix, [_label(ba) for ba in ba_sizes], [_label(ab) for ab in ab_sizes]


def _ab_ba_family_sizes(module: BaseMultiqcModule) -> Set[str]:
    rows_by_sample: Dict[str, List[DuplexFamilySizeMetric]] = {}
    for f in module.find_log_files("fgumi/duplex_family_sizes"):
        rows = load_rows(f, DuplexFamilySizeMetric)
        if rows is None:
            continue
        s_name = sample_name(module, f)
        module.add_data_source(f, s_name)
        module.add_software_version(None, s_name)
        rows_by_sample[s_name] = rows
    rows_by_sample = module.ignore_samples(rows_by_sample)
    if not rows_by_sample:
        return set()

    min_strand: Dict[str, Dict[int, int]] = {}
    for s_name, rows in rows_by_sample.items():
        totals: Dict[int, int] = defaultdict(int)
        for r in rows:
            totals[min(r.ab_size, r.ba_size)] += r.count
        min_strand[s_name] = dict(sorted(totals.items()))
    module.add_section(
        name="Duplex family sizes: smaller strand",
        anchor="fgumi-duplex-min-strand",
        description="Duplex families by the size of their smaller strand, min(AB, BA), across samples.",
        helptext="A family with a smaller strand of size 0 has reads from only one strand and cannot form a duplex.",
        plot=linegraph.plot(
            min_strand,
            {
                "id": "fgumi_duplex_min_strand",
                "title": "fgumi: Duplex families by smaller strand",
                "xlab": "min(AB reads, BA reads)",
                "ylab": "Families",
            },
        ),
    )
    for s_name, rows in rows_by_sample.items():
        for suffix, min_ba, title in [("", 0, "all families"), ("_both", 1, "families with AB > 0 and BA > 0")]:
            built = _heatmap_matrix(rows, min_ba)
            if built is None:
                continue
            matrix, xcats, ycats = built
            plot_id = f"fgumi_duplex_ab_ba{suffix}_{safe_id(s_name)}"
            module.add_section(
                name=f"Duplex family sizes: {s_name} ({title})",
                anchor=plot_id.replace("_", "-"),
                description="log10(families) by AB and BA strand size; sizes of 20 or more are pooled as 20+.",
                plot=heatmap.plot(
                    matrix,
                    xcats,
                    ycats,
                    {
                        "id": plot_id,
                        "title": f"fgumi: Duplex AB x BA family sizes, {s_name}",
                        "xlab": "BA reads",
                        "ylab": "AB reads",
                        "zlab": "log10(families)",
                        "xcats_samples": False,
                        "ycats_samples": False,
                        "cluster_rows": False,
                        "cluster_cols": False,
                    },
                ),
            )
    module.write_data_file(min_strand, "multiqc_fgumi_duplex_min_strand")
    return set(rows_by_sample)
