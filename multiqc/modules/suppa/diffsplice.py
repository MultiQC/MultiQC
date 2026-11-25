import logging
import math
from collections.abc import Iterable

from multiqc import BaseMultiqcModule, config
from multiqc.plots import bargraph, table
from multiqc.plots.bargraph import BarPlotConfig
from multiqc.plots.table_object import ColumnDict, TableConfig
from multiqc.types import SectionAlert

from .events import EVENT_TYPE_COLOURS, ISOFORM, event_type

log = logging.getLogger(__name__)

# The defaults of SUPPA `clusterEvents`, which picks significant events with
# p < --sig-threshold and |dPSI| > --dpsi-threshold
DEFAULT_PVALUE_CUTOFF = 0.05
DEFAULT_DPSI_CUTOFF = 0.05


def parse_dpsi(lines: Iterable[str], pvalue_cutoff: float, dpsi_cutoff: float) -> dict[str, dict]:
    """
    Parse a SUPPA `diffSplice` `.dpsi` file.

    The header holds a `<cond1>-<cond2>_dPSI` and `<cond1>-<cond2>_p-val` column pair for
    each comparison, and each row an event id followed by those values. Returns per
    comparison counts of events, tested events and significant events, overall and by
    event type.
    """
    line_iter = iter(lines)
    try:
        header = next(line_iter).rstrip("\r\n").split("\t")
    except StopIteration:
        raise ValueError("file is empty")
    if not header or len(header) % 2 != 0:
        raise ValueError(f"expected dPSI and p-value column pairs, got {len(header)} columns")
    comparisons = []
    for dpsi_col, pval_col in zip(header[::2], header[1::2]):
        if not dpsi_col.endswith("_dPSI") or pval_col != dpsi_col[: -len("_dPSI")] + "_p-val":
            raise ValueError(f"unexpected column pair '{dpsi_col}', '{pval_col}'")
        comparisons.append(dpsi_col[: -len("_dPSI")])

    data: dict[str, dict] = {
        c: {
            "events": 0,
            "tested": 0,
            "significant": 0,
            "significant_up": 0,
            "significant_down": 0,
            "tested_by_type": {},
            "significant_by_type": {},
        }
        for c in comparisons
    }
    for line_no, line in enumerate(line_iter, start=2):
        fields = line.rstrip("\r\n").split("\t")
        if fields == [""]:
            continue
        if len(fields) != len(header) + 1:
            raise ValueError(f"line {line_no} has {len(fields) - 1} values, expected {len(header)}")
        etype = event_type(fields[0])
        for i, comparison in enumerate(comparisons):
            d = data[comparison]
            dpsi = float(fields[1 + 2 * i])
            pval = float(fields[2 + 2 * i])
            d["events"] += 1
            # Events that were not tested, or absent from this comparison, have dPSI nan
            if math.isnan(dpsi) or math.isnan(pval):
                continue
            d["tested"] += 1
            d["tested_by_type"][etype] = d["tested_by_type"].get(etype, 0) + 1
            if pval < pvalue_cutoff and abs(dpsi) > dpsi_cutoff:
                d["significant"] += 1
                d["significant_up" if dpsi > 0 else "significant_down"] += 1
                d["significant_by_type"][etype] = d["significant_by_type"].get(etype, 0) + 1

    for d in data.values():
        if d["events"]:
            d["pct_tested"] = 100.0 * d["tested"] / d["events"]
        if d["tested"]:
            d["pct_significant"] = 100.0 * d["significant"] / d["tested"]
    return data


def parse_suppa_diffsplice(module: BaseMultiqcModule) -> int:
    """Find SUPPA `.dpsi` files and add differential splicing sections. Returns the number of comparisons."""
    cfg = getattr(config, "suppa", {})
    pvalue_cutoff = float(cfg.get("pvalue_cutoff", DEFAULT_PVALUE_CUTOFF))
    dpsi_cutoff = float(cfg.get("dpsi_cutoff", DEFAULT_DPSI_CUTOFF))

    data: dict[str, dict] = {}
    for f in module.find_log_files("suppa/diffsplice", filehandles=True):
        try:
            parsed = parse_dpsi(f["f"], pvalue_cutoff, dpsi_cutoff)
        except ValueError as e:
            log.warning(f"Skipping {f['fn']}: {e}")
            continue
        for comparison, comparison_data in parsed.items():
            # The output prefix names the run. Files with several comparisons get the
            # comparison label too.
            s_name = f["s_name"] if len(parsed) == 1 else f"{f['s_name']} {comparison}"
            if s_name in data:
                log.debug(f"Duplicate comparison name found in {f['fn']}! Overwriting: {s_name}")
            comparison_data["comparison"] = comparison
            module.add_data_source(f, s_name=s_name, section="diffsplice")
            module.add_software_version(None, s_name)
            data[s_name] = comparison_data

    data = module.ignore_samples(data)
    if not data:
        return 0

    _add_summary_section(module, data, pvalue_cutoff, dpsi_cutoff)
    _add_significant_by_type_section(module, data, pvalue_cutoff, dpsi_cutoff)

    module.write_data_file(
        {
            s: {
                **{k: v for k, v in d.items() if not k.endswith("_by_type")},
                **{f"tested {t}": n for t, n in d["tested_by_type"].items()},
                **{f"significant {t}": n for t, n in d["significant_by_type"].items()},
            }
            for s, d in data.items()
        },
        "multiqc_suppa_diffsplice",
    )
    return len(data)


def _significance_text(pvalue_cutoff: float, dpsi_cutoff: float) -> str:
    return f"p < {pvalue_cutoff:g} and |dPSI| > {dpsi_cutoff:g}"


def _add_summary_section(
    module: BaseMultiqcModule, data: dict[str, dict], pvalue_cutoff: float, dpsi_cutoff: float
) -> None:
    significance = _significance_text(pvalue_cutoff, dpsi_cutoff)
    headers: dict[str, ColumnDict] = {
        "comparison": {
            "title": "Conditions",
            "description": "Conditions compared, from the dpsi file header: `<condition 1>-<condition 2>`",
            "scale": False,
        },
        "events": {
            "title": "Events",
            "description": "Events in the dpsi file",
            "format": "{:,.0f}",
            "scale": "Greys",
            "hidden": True,
        },
        "tested": {
            "title": "Tested",
            "description": "Events with a dPSI value. Events with missing PSI values in a condition (see `--nan-threshold`) are not tested.",
            "format": "{:,.0f}",
            "scale": "Blues",
        },
        "pct_tested": {
            "title": "% Tested",
            "description": "Percentage of events with a dPSI value",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "GnBu",
            "hidden": True,
        },
        "significant": {
            "title": "Significant",
            "description": f"Tested events with {significance}",
            "format": "{:,.0f}",
            "scale": "Purples",
        },
        "pct_significant": {
            "title": "% Significant",
            "description": f"Percentage of the tested events that are significant ({significance})",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdPu",
        },
        "significant_up": {
            "title": "Higher PSI",
            "description": f"Significant events ({significance}) with dPSI > 0, a higher PSI in condition 2 "
            "than in condition 1",
            "format": "{:,.0f}",
            "scale": "Reds",
        },
        "significant_down": {
            "title": "Lower PSI",
            "description": f"Significant events ({significance}) with dPSI < 0, a lower PSI in condition 2 "
            "than in condition 1",
            "format": "{:,.0f}",
            "scale": "Blues",
        },
    }
    module.add_section(
        name="Differential Splicing",
        anchor="suppa_differential_splicing",
        description=f"Events tested and called significant in each SUPPA `diffSplice` comparison. "
        f"Significant: {significance}.",
        helptext=f"""
Each row is one comparison from a `.dpsi` file, named after the file (the `-o` output
prefix). A file with several comparisons, from three or more conditions, gives one row
per comparison.

Direction: SUPPA `dPSI` is the PSI of condition 2 minus that of condition 1 (the second
and first conditions in the label), so **Higher PSI** means `dPSI` > 0.

Events are counted as significant with {significance}, the defaults of SUPPA
`clusterEvents`. The p-values are those written by SUPPA, so they are corrected for
multiple testing only if `diffSplice` ran with `--gene-correction`. Set the thresholds in
a MultiQC config file (defaults shown):

```yaml
suppa:
  pvalue_cutoff: {DEFAULT_PVALUE_CUTOFF:g}
  dpsi_cutoff: {DEFAULT_DPSI_CUTOFF:g}
```
""",
        plot=table.plot(
            data,
            headers,
            TableConfig(
                id="suppa_differential_splicing_table",
                title="SUPPA: Differential Splicing",
                col1_header="Comparison",
            ),
        ),
    )


def _add_significant_by_type_section(
    module: BaseMultiqcModule, data: dict[str, dict], pvalue_cutoff: float, dpsi_cutoff: float
) -> None:
    # Isoform comparisons have a single category, which the summary table already shows
    local = {s: d for s, d in data.items() if set(d["tested_by_type"]) - {ISOFORM}}
    if not local:
        return
    cats = {t: {"name": t, "color": colour} for t, colour in EVENT_TYPE_COLOURS.items()}

    no_significant = sorted(s for s, d in local.items() if d["significant"] == 0)
    alerts = []
    plot = None
    if len(no_significant) == len(local):
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
            {s: d["significant_by_type"] for s, d in local.items()},
            cats,
            BarPlotConfig(
                id="suppa_significant_events_plot",
                title="SUPPA: Significant Events",
                ylab="Significant events",
                cpswitch_counts_label="Number of events",
            ),
        )

    module.add_section(
        name="Significant Events",
        anchor="suppa_significant_events",
        description="Significant events in each comparison of local events (`psiPerEvent`), by event type. "
        f"Significant: {_significance_text(pvalue_cutoff, dpsi_cutoff)}.",
        helptext="""
SUPPA event types:

* **SE**: skipped exon
* **A5** / **A3**: alternative 5' / 3' splice site
* **MX**: mutually exclusive exons
* **RI**: retained intron
* **AF** / **AL**: alternative first / last exon

Comparisons of isoforms (`psiPerIsoform`) are not shown, as they have no event types.
""",
        alerts=alerts,
        plot=plot,
    )
