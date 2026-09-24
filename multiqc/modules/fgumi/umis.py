"""Per-UMI files, summarized rather than listed: `<prefix>.umi_counts.txt` (simplex/duplex),
`<prefix>.duplex_umi_counts.txt`, and `correct --metrics`. Rows are streamed so a large file is not turned into
one pydantic object per UMI."""

import itertools
from collections import Counter
from typing import Any, Dict, Iterable, Set

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph

from .schemas import DuplexUmiMetric, MetricFormatError, UmiCorrectionMetric, UmiMetric
from .util import drop_none, register, skip_unreadable, stream_dicts


def summarize_observations(counts: Iterable[int]) -> Dict[str, Any]:
    """Number of UMIs, median raw observations, % seen once, and UMIs per log2 bin of raw observations.

    Tallies the distinct observation counts rather than keeping one value per UMI, so memory stays small for a
    file with millions of UMIs.
    """
    tally: Counter = Counter(counts)
    n = sum(tally.values())
    bins: Counter = Counter()
    for value, umis in tally.items():
        bins[1 << (value - 1).bit_length() if value > 1 else 1] += umis
    # With no UMIs the median and singleton share are undefined, not zero.
    return {
        "n": n,
        "median": _median(tally, n) if n else None,
        "singleton_pct": 100.0 * tally[1] / n if n else None,
        "bins": dict(sorted(bins.items())),
    }


def _median(tally: Counter, n: int) -> float:
    """The median of ``n`` values given as ``{value: occurrences}``."""
    ordered = sorted(tally.items())

    def nth(index: int) -> int:
        seen = 0
        for value, occurrences in ordered:
            seen += occurrences
            if index < seen:
                return value
        raise IndexError(index)

    return float(nth(n // 2)) if n % 2 else (nth(n // 2 - 1) + nth(n // 2)) / 2.0


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
            schema = DuplexUmiMetric if DuplexUmiMetric.matches(header_line.rstrip("\r\n").split("\t")) else UmiMetric
            kind = "duplex_umi_counts" if schema is DuplexUmiMetric else "umi_counts"
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
                {s: data[s]["bins"] for s in data},
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
                "scale": "Purples",
                "hidden": True,
                "description": f"Distinct UMIs ({title.lower()})",
                "format": "{:,.0f}",
            },
            "median": {
                "title": "Median obs",
                "scale": "Blues",
                "hidden": True,
                "description": "Median raw observations per UMI",
            },
            "singleton_pct": {
                "title": "% singleton",
                "scale": "YlOrRd",
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
