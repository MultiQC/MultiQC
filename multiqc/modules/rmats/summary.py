import logging
import re

from multiqc import BaseMultiqcModule
from multiqc.plots import bargraph, table
from multiqc.plots.table_object import ColumnDict
from multiqc.types import SectionAlert

log = logging.getLogger(__name__)

EVENT_TYPES: dict[str, dict[str, str]] = {
    "SE": {"name": "Skipped exon", "color": "#437bb1", "light": "#a6c4e2"},
    "A5SS": {"name": "Alternative 5' splice site", "color": "#e8853a", "light": "#f5c9a4"},
    "A3SS": {"name": "Alternative 3' splice site", "color": "#3a9a5b", "light": "#a8d8b5"},
    "MXE": {"name": "Mutually exclusive exons", "color": "#8e5ea2", "light": "#d1b8dd"},
    "RI": {"name": "Retained intron", "color": "#c0392b", "light": "#eba9a2"},
}

COUNT_COLUMNS = [
    "TotalEventsJC",
    "TotalEventsJCEC",
    "SignificantEventsJC",
    "SigEventsJCSample1HigherInclusion",
    "SigEventsJCSample2HigherInclusion",
    "SignificantEventsJCEC",
    "SigEventsJCECSample1HigherInclusion",
    "SigEventsJCECSample2HigherInclusion",
]

# Added to the file name by `rMATS_P/summary.py --summary-prefix`
CUTOFF_SUFFIX = re.compile(r"_(FDR|PValue)_([^_]+)_IncLevelDifference_([^_]+)\.txt$")

COUNTS = ["JC", "JCEC"]

# The cutoff rmats.py passes to rMATS_P/summary.py when it writes summary.txt
DEFAULT_SIGNIFICANCE = "FDR <= 0.05 (rMATS default)"


def parse_summary_file(lines) -> dict[str, int]:
    """
    Parse an rMATS-turbo summary file into {"<event type>_<column>": count}.
    The EventTypeDescription column, added in v4.2.0, is ignored.
    """
    header: list[str] | None = None
    data: dict[str, int] = {}
    for line in lines:
        line = line.rstrip("\r\n")
        if not line.strip():
            continue
        cols = line.split("\t")
        if header is None:
            header = cols
            missing = [c for c in ["EventType"] + COUNT_COLUMNS if c not in header]
            if missing:
                raise ValueError(f"summary header is missing columns: {', '.join(missing)}")
            continue
        if len(cols) != len(header):
            raise ValueError(f"expected {len(header)} columns, found {len(cols)} in line '{line}'")
        row = dict(zip(header, cols))
        event = row["EventType"]
        if event not in EVENT_TYPES:
            raise ValueError(f"unknown event type '{event}'")
        for col in COUNT_COLUMNS:
            data[f"{event}_{col}"] = int(row[col])

    missing_events = [e for e in EVENT_TYPES if f"{e}_TotalEventsJC" not in data]
    if missing_events:
        raise ValueError(f"summary has no rows for event types: {', '.join(missing_events)}")
    return data


def summary_name_and_cutoff(fn: str, root: str) -> tuple[str, str | None]:
    """
    Comparison name, before sample name cleaning, and significance cutoff from the summary
    file name. A cutoff suffix stays in the name, so that summaries of one comparison at
    several cutoffs remain distinct. A file named only `summary` takes the name of its
    directory.
    """
    cutoff = None
    suffix = ""
    base = fn
    m = CUTOFF_SUFFIX.search(fn)
    if m:
        cutoff = f"{m.group(1)} <= {m.group(2)}, |IncLevelDifference| >= {m.group(3)}"
        suffix = fn[m.start() : -len(".txt")]
        base = fn[: m.start()]
    if re.fullmatch(r"summary(\.txt)?", base, flags=re.IGNORECASE):
        dirs = [d for d in re.split(r"[\\/]", root) if d and d != "."]
        base = dirs[-1] if dirs else "rmats"
    return base + suffix, cutoff


def parse_rmats_summary(module: BaseMultiqcModule) -> int:
    data: dict[str, dict[str, int]] = {}
    cutoffs: dict[str, str] = {}

    for f in module.find_log_files("rmats/summary", filehandles=True):
        try:
            parsed = parse_summary_file(f["f"])
        except ValueError as e:
            log.warning(f"Skipping {f['fn']}: {e}")
            continue
        name, cutoff = summary_name_and_cutoff(f["fn"], f["root"])
        s_name = module.clean_s_name(name, f)
        if s_name in data:
            log.debug(f"Duplicate comparison name found in {f['fn']}! Overwriting: {s_name}")
        data[s_name] = parsed
        if cutoff is not None:
            cutoffs[s_name] = cutoff
        module.add_data_source(f, s_name=s_name, section="summary")
        module.add_software_version(None, s_name)

    data = module.ignore_samples(data)
    if len(data) == 0:
        return 0

    _add_counts_table(module, data, cutoffs)
    _add_significant_section(module, data)

    write_data: dict[str, dict[str, int | str]] = {s: dict(d) for s, d in data.items()}
    for s_name, cutoff in cutoffs.items():
        if s_name in write_data:
            write_data[s_name]["cutoff"] = cutoff
    module.write_data_file(write_data, "multiqc_rmats_summary")

    return len(data)


def _add_counts_table(module: BaseMultiqcModule, data: dict[str, dict[str, int]], cutoffs: dict[str, str]):
    rows: dict[str, dict[str, int | float | str]] = {}
    for s_name, d in data.items():
        row: dict[str, int | float | str] = {}
        for counts in COUNTS:
            tested = sum(d[f"{e}_TotalEvents{counts}"] for e in EVENT_TYPES)
            significant = sum(d[f"{e}_SignificantEvents{counts}"] for e in EVENT_TYPES)
            row[f"tested_{counts}"] = tested
            row[f"significant_{counts}"] = significant
            if tested:
                row[f"significant_pct_{counts}"] = significant / tested * 100.0
        if s_name in cutoffs:
            row["cutoff"] = cutoffs[s_name]
        rows[s_name] = row

    headers: dict[str, ColumnDict] = {}
    for counts in COUNTS:
        headers[f"tested_{counts}"] = {
            "title": f"Tested ({counts})",
            "description": f"Events tested for differential splicing with {counts} counts",
            "format": "{:,.0f}",
            "scale": "Blues",
        }
        headers[f"significant_{counts}"] = {
            "title": f"Significant ({counts})",
            "description": f"Events called significant with {counts} counts",
            "format": "{:,.0f}",
            "scale": "Purples",
        }
        headers[f"significant_pct_{counts}"] = {
            "title": f"% Significant ({counts})",
            "description": f"Percentage of the events tested with {counts} counts that are significant",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdPu",
        }
    if cutoffs:
        headers["cutoff"] = {
            "title": "Cutoff",
            "description": "Significance cutoff, from the summary file name written by rMATS_P/summary.py "
            "--summary-prefix",
            "scale": False,
        }

    module.add_section(
        name="Event Counts",
        anchor="rmats_event_counts",
        description=f"Splicing events tested and called significant in each rMATS comparison, over all event "
        f"types. Significant: {DEFAULT_SIGNIFICANCE}, unless the Cutoff column gives other cutoffs.",
        helptext="""
        Each row is one rMATS comparison of the `--b1` group against the `--b2` group.
        rMATS counts each event in two ways: **JC** uses only reads that span splice junctions,
        **JCEC** also uses reads on the exons themselves. An event is tested when each group has
        reads for it and the samples have both inclusion and skipping reads. Few tested events
        usually means low depth, or reads that rMATS discarded (see the Read Outcomes section).

        The significance cutoffs are set when rMATS writes the summary, not by MultiQC.
        `rMATS_P/summary.py` can recount with other cutoffs; when it writes the cutoffs into the
        file name they are shown in the Cutoff column.
        """,
        plot=table.plot(
            rows,
            headers,
            {
                "id": "rmats_event_counts_table",
                "title": "rMATS: Event Counts",
                "namespace": "rMATS",
                "col1_header": "Comparison",
            },
        ),
    )


def _add_significant_section(module: BaseMultiqcModule, data: dict[str, dict[str, int]]):
    cats: dict[str, dict[str, str]] = {}
    for event, info in EVENT_TYPES.items():
        cats[f"{event}_1"] = {"name": f"{info['name']}: higher inclusion in group 1", "color": info["color"]}
        cats[f"{event}_2"] = {"name": f"{info['name']}: higher inclusion in group 2", "color": info["light"]}

    datasets = []
    for counts in COUNTS:
        datasets.append(
            {
                s_name: {
                    f"{event}_{group}": d[f"{event}_SigEvents{counts}Sample{group}HigherInclusion"]
                    for event in EVENT_TYPES
                    for group in (1, 2)
                }
                for s_name, d in data.items()
            }
        )

    no_significant = sorted(
        s for s, d in data.items() if all(d[f"{e}_SignificantEvents{c}"] == 0 for e in EVENT_TYPES for c in COUNTS)
    )
    alerts = []
    plot = None
    if len(no_significant) == len(data):
        alerts.append(SectionAlert(message="No comparison has any significant events.", level="info"))
    else:
        if no_significant:
            alerts.append(
                SectionAlert(
                    message="These comparisons have no significant events.",
                    level="info",
                    affected_samples=no_significant,
                )
            )
        plot = bargraph.plot(
            datasets,
            [cats, cats],
            {
                "id": "rmats_significant_events_plot",
                "title": "rMATS: Significant Events",
                "ylab": "Significant events",
                "cpswitch": False,
                "data_labels": [
                    {"name": "JC", "ylab": "Significant events (JC)"},
                    {"name": "JCEC", "ylab": "Significant events (JCEC)"},
                ],
            },
        )

    module.add_section(
        name="Significant Events",
        anchor="rmats_significant_events",
        description="Significant differential splicing events in each comparison, by event type and by the "
        f"group with higher inclusion. Significant: {DEFAULT_SIGNIFICANCE}, unless the summary file name "
        "gives other cutoffs.",
        helptext="""
        Group 1 is the `--b1` (or `--s1`) group and group 2 the `--b2` (or `--s2`) group.
        Direction: rMATS `IncLevelDifference` is the inclusion level of group 1 minus that of
        group 2, so **higher inclusion in group 1** means `IncLevelDifference` > 0. For a skipped
        exon, the exon is included more often in group 1.

        * **SE**: skipped exon
        * **A5SS**: alternative 5' splice site
        * **A3SS**: alternative 3' splice site
        * **MXE**: mutually exclusive exons
        * **RI**: retained intron

        Use the buttons to switch between JC counts (junction reads only) and JCEC counts
        (junction and exon reads).
        """,
        alerts=alerts,
        plot=plot,
    )
