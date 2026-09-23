"""Per-UMI files, summarized rather than listed: `<prefix>.umi_counts.txt` (simplex/duplex),
`<prefix>.duplex_umi_counts.txt`, and `correct --metrics`. Rows are streamed so a large file is not turned into
one pydantic object per UMI."""

import itertools
import statistics
from collections import Counter
from typing import Any, Dict, Iterable, List, Set

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph

from .schemas import DuplexUmiMetric, MetricFormatError, UmiCorrectionMetric, UmiMetric
from .util import drop_none, register, skip_unreadable, stream_dicts


def summarize_observations(counts: Iterable[int]) -> Dict[str, Any]:
    """Number of UMIs, median raw observations, % seen once, and UMIs per log2 bin of raw observations."""
    values: List[int] = list(counts)
    bins: Counter = Counter(1 << (value - 1).bit_length() if value > 1 else 1 for value in values)
    # With no UMIs the median and singleton share are undefined, not zero.
    return {
        "n": len(values),
        "median": float(statistics.median(values)) if values else None,
        "singleton_pct": 100.0 * sum(1 for v in values if v == 1) / len(values) if values else None,
        "bins": dict(sorted(bins.items())),
    }


def parse_reports(module: BaseMultiqcModule) -> Set[str]:
    return _umi_counts(module) | _correct(module)


def _umi_counts(module: BaseMultiqcModule) -> Set[str]:
    summaries: Dict[str, Dict[str, Dict[str, Any]]] = {"umi_counts": {}, "duplex_umi_counts": {}}
    for f in module.find_log_files("fgumi/umi_counts", filehandles=True):
        handle: Any = f["f"]  # a text handle: found with filehandles=True
        if handle is None:
            continue
        try:
            header_line = handle.readline()
            kind = "duplex_umi_counts" if "fraction_unique_observations_expected" in header_line else "umi_counts"
            schema = DuplexUmiMetric if kind == "duplex_umi_counts" else UmiMetric
            rows = stream_dicts(itertools.chain([header_line], handle), f["fn"], schema)
            summary = summarize_observations(int(row["raw_observations"]) for row in rows)
        except (MetricFormatError, ValueError) as error:  # UnicodeDecodeError is a ValueError
            skip_unreadable(f, error)
            continue
        summaries[kind][register(module, f, section=kind)] = summary

    parsed: Set[str] = set()
    for kind, title in [("umi_counts", "UMI counts"), ("duplex_umi_counts", "Duplex UMI counts")]:
        data = module.ignore_samples(summaries[kind])
        if not data:
            continue
        module.add_section(
            name=title,
            anchor=f"fgumi-{kind.replace('_', '-')}",
            description="How many UMIs were seen, binned by how many read pairs carried each one (log2 bins).",
            helptext="Summarized from the per-UMI file written by `fgumi simplex-metrics` / `duplex-metrics` "
            "(or fgbio CollectDuplexSeqMetrics); "
            "individual UMIs are not listed.",
            plot=linegraph.plot(
                {s: dict(data[s]["bins"]) for s in data},
                {
                    "id": f"fgumi_{kind}",
                    "title": f"fgumi: {title}",
                    "xlab": "Raw observations per UMI (log2 bin upper bound)",
                    "ylab": "UMIs",
                    "xlog": True,
                },
            ),
        )
        table = {s: drop_none({k: v for k, v in data[s].items() if k != "bins"}) for s in data}
        headers: Dict[str, Dict[str, Any]] = {
            "n": {
                "title": "UMIs" if kind == "umi_counts" else "Duplex UMIs",
                "hidden": True,
                "description": f"Distinct UMIs ({title.lower()})",
                "format": "{:,.0f}",
            },
            "median": {"title": "Median obs", "hidden": True, "description": "Median raw observations per UMI"},
            "singleton_pct": {
                "title": "% singleton",
                "hidden": True,
                "suffix": "%",
                "max": 100,
                "min": 0,
                "description": "UMIs seen in exactly one read pair (raw observation)",
            },
        }
        module.general_stats_addcols(table, headers, namespace=title)
        module.write_data_file(table, f"multiqc_fgumi_{kind}")
        parsed |= set(data)
    return parsed


def _correct(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[str, int]] = {}
    for f in module.find_log_files("fgumi/correct", filehandles=True):
        if f["f"] is None:
            continue
        totals = {"perfect": 0, "one_mismatch": 0, "two_mismatch": 0, "other": 0, "unmatched": 0}
        try:
            for row in stream_dicts(f["f"], f["fn"], UmiCorrectionMetric):
                if set(row["umi"].replace("-", "")) == {"N"}:
                    totals["unmatched"] += int(row["total_matches"])
                    continue
                totals["perfect"] += int(row["perfect_matches"])
                totals["one_mismatch"] += int(row["one_mismatch_matches"])
                totals["two_mismatch"] += int(row["two_mismatch_matches"])
                totals["other"] += int(row["other_matches"])
        except (MetricFormatError, ValueError) as error:  # UnicodeDecodeError is a ValueError
            skip_unreadable(f, error)
            continue
        data[register(module, f)] = totals
    data = module.ignore_samples(data)
    if not data:
        return set()
    module.add_section(
        name="UMI correction",
        anchor="fgumi-correct",
        description="UMI observations that matched an expected UMI exactly, after correcting mismatches, or not at "
        "all.",
        helptext="Written by `fgumi correct --metrics` (or fgbio CorrectUmis). Each UMI segment of each read counts "
        "once, so a dual UMI counts twice; the all-N row counts segments that matched no expected UMI.",
        plot=bargraph.plot(
            data,
            {
                "perfect": {"name": "Perfect match"},
                "one_mismatch": {"name": "1 mismatch"},
                "two_mismatch": {"name": "2 mismatches"},
                "other": {"name": "3+ mismatches"},
                "unmatched": {"name": "Unmatched"},
            },
            {"id": "fgumi_correct", "title": "fgumi: UMI correction", "ylab": "UMI observations"},
        ),
    )
    module.write_data_file(data, "multiqc_fgumi_correct")
    return set(data)
