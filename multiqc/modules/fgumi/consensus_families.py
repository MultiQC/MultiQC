"""Consensus family sizes from `simplex-metrics` / `duplex-metrics` (and the consensus callers' `--metrics`):
per-strand family sizes (`<prefix>.family_sizes.txt`) and the duplex AB x BA family sizes
(`<prefix>.duplex_family_sizes.txt`), drawn like fgbio CollectDuplexSeqMetrics' plots #1-#3."""

import html
import logging
import math
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

from multiqc import config
from multiqc.base_module import BaseMultiqcModule
from multiqc.report import clean_htmlid
from multiqc.plots import heatmap, linegraph
from multiqc.types import SectionAlert

from .schemas import DuplexFamilySizeMetric, DuplexStrandFamilySizeMetric, SimplexFamilySizeMetric
from .util import drop_none, flatten, header_columns, iter_samples, load_rows, pct, register

log = logging.getLogger(__name__)

HEATMAP_CAP = 20  # sizes at or above this fold into one "20+" row/column
# Above this many samples the per-sample heatmaps are skipped (they would swamp the report) and only the
# cross-sample min(AB, BA) line is drawn. Override with `fgumi_config: {max_duplex_heatmap_samples: N}`.
MAX_HEATMAP_SAMPLES = 5

_HEATMAP_HELP = (
    "Each cell counts duplex families by the number of reads from the AB and BA strands, as in fgbio "
    "CollectDuplexSeqMetrics. Families on the diagonal are balanced; a family with 0 reads on one strand cannot "
    "form a duplex."
)


def parse_reports(module: BaseMultiqcModule) -> Set[str]:
    return _strand_family_sizes(module) | _ab_ba_family_sizes(module)


def _strand_family_sizes(module: BaseMultiqcModule) -> Set[str]:
    by_kind: Dict[str, Dict[str, Dict[str, Dict[int, Optional[float]]]]] = {"simplex": {}, "duplex": {}}
    for f in module.find_log_files("fgumi/strand_family_sizes"):
        if f["f"] is None:
            continue
        is_duplex = DuplexStrandFamilySizeMetric.matches(header_columns(f))
        rows = load_rows(f, DuplexStrandFamilySizeMetric if is_duplex else SimplexFamilySizeMetric)
        if rows is None:
            continue
        # Not "<kind>_family_sizes": "duplex_family_sizes" is already the section of the AB x BA files.
        s_name = register(module, f, section=f"strand_family_sizes_{'duplex' if is_duplex else 'simplex'}")
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
            description=f"Family size distributions from `fgumi {kind}-metrics`"
            + (" (or fgbio CollectDuplexSeqMetrics)" if kind == "duplex" else "")
            + ", one tab per grouping.",
            helptext="CS: grouped by coordinate and strand only. SS: also by UMI (single-strand families)."
            + (
                " DS: single-strand families paired into double-strand families. Matches fgbio CollectDuplexSeqMetrics."
                if kind == "duplex"
                else ""
            ),
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
            flatten({s: {strand: data[s][strand] for strand in strands} for s in data}),
            f"multiqc_fgumi_{kind}_family_sizes",
        )
        parsed |= set(data)
    return parsed


def _label(size: int) -> str:
    return f"{HEATMAP_CAP}+" if size >= HEATMAP_CAP else str(size)


def _heatmap_matrix(
    rows: List[DuplexFamilySizeMetric], min_ba: int
) -> Optional[Tuple[List[List[Optional[float]]], List[str], List[str]]]:
    """log10(count) with a row per capped BA size and a column per capped AB size (fgbio's x = AB, y = BA), and
    the column and row labels; rows with ``ba_size < min_ba`` are dropped."""
    cells: Dict[Tuple[int, int], int] = defaultdict(int)
    for r in rows:
        if r.ba_size >= min_ba:
            cells[(min(r.ab_size, HEATMAP_CAP), min(r.ba_size, HEATMAP_CAP))] += r.count
    if not cells:
        return None
    # Contiguous axes, so a size with no families is an empty cell rather than missing from the axis.
    ab_sizes = list(range(min(ab for ab, _ in cells), max(ab for ab, _ in cells) + 1))
    ba_sizes = list(range(min(ba for _, ba in cells), max(ba for _, ba in cells) + 1))
    matrix = [[math.log10(cells[(ab, ba)]) if cells.get((ab, ba)) else None for ab in ab_sizes] for ba in ba_sizes]
    return matrix, [_label(ab) for ab in ab_sizes], [_label(ba) for ba in ba_sizes]


def _ab_ba_family_sizes(module: BaseMultiqcModule) -> Set[str]:
    rows_by_sample = dict(iter_samples(module, "fgumi/duplex_family_sizes", DuplexFamilySizeMetric))
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
    max_samples = (getattr(config, "fgumi_config", None) or {}).get("max_duplex_heatmap_samples", MAX_HEATMAP_SAMPLES)
    if len(rows_by_sample) > max_samples:
        log.info(f"Skipping per-sample duplex heatmaps for {len(rows_by_sample)} samples (limit {max_samples})")
    else:
        _per_sample_heatmaps(module, rows_by_sample)
    module.write_data_file(min_strand, "multiqc_fgumi_duplex_min_strand")
    return set(rows_by_sample)


def _per_sample_heatmaps(module: BaseMultiqcModule, rows_by_sample: Dict[str, List[DuplexFamilySizeMetric]]) -> None:
    for s_name, rows in rows_by_sample.items():
        if not rows:
            module.add_section(
                name=f"Duplex family sizes: {html.escape(s_name)}",
                anchor=clean_htmlid(f"fgumi_duplex_ab_ba_{s_name}").replace("_", "-"),
                description="No duplex families were recorded.",
                helptext=_HEATMAP_HELP,
                plot=None,
                alerts=[SectionAlert(message="No duplex families were recorded", affected_samples=[s_name])],
            )
            continue
        variants = [("", 0, "all families"), ("_both", 1, "families with AB > 0 and BA > 0")]
        if all(r.ba_size > 0 for r in rows):
            variants = variants[:1]  # every family has BA reads, so the both-strands heatmap would be identical
        for suffix, min_ba, title in variants:
            built = _heatmap_matrix(rows, min_ba)
            plot_id = clean_htmlid(f"fgumi_duplex_ab_ba{suffix}_{s_name}")
            if built is None:
                # "No families with reads on both strands" is exactly what this plot exists to show.
                module.add_section(
                    name=f"Duplex family sizes: {html.escape(s_name)} ({title})",
                    anchor=plot_id.replace("_", "-"),
                    description="No families had reads from both strands, so no duplexes can be formed.",
                    helptext=_HEATMAP_HELP,
                    plot=None,
                    alerts=[SectionAlert(message="No families with reads on both strands", affected_samples=[s_name])],
                )
                continue
            matrix, xcats, ycats = built
            module.add_section(
                name=f"Duplex family sizes: {html.escape(s_name)} ({title})",
                anchor=plot_id.replace("_", "-"),
                description=f"log10(families) by AB and BA strand size; sizes of {HEATMAP_CAP} or more are pooled "
                f"as {HEATMAP_CAP}+.",
                helptext=_HEATMAP_HELP,
                plot=heatmap.plot(
                    matrix,
                    xcats,
                    ycats,
                    {
                        "id": plot_id,
                        "title": f"fgumi: Duplex AB x BA family sizes, {html.escape(s_name)}",
                        "xlab": "AB reads",
                        "ylab": "BA reads",
                        "zlab": "log10(families)",
                        "xcats_samples": False,
                        "ycats_samples": False,
                        "cluster_rows": False,
                        "cluster_cols": False,
                    },
                ),
            )
