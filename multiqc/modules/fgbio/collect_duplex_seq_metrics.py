import logging
from typing import Dict, List

from multiqc import BaseMultiqcModule
from multiqc.plots import linegraph

log = logging.getLogger(__name__)

MAX_FAMILY_SIZE = 20


def run_collect_duplex_seq_metrics(module: BaseMultiqcModule) -> int:
    """Parse output from fgbio CollectDuplexSeqMetrics (family_sizes.txt)."""
    data = _parse_family_sizes(module)
    if not data:
        return 0

    _add_family_size_plot(module, data)
    _add_general_stats(module, data)
    return len(data)


def _parse_family_sizes(module: BaseMultiqcModule) -> Dict[str, List[Dict]]:
    """Parse family_sizes.txt files produced by CollectDuplexSeqMetrics.

    Each row has columns: family_size, cs_count, cs_fraction, ...,
    ss_count, ss_fraction, ..., ds_count, ds_fraction, ...

    Orphaned SSCS count = ss_count - ds_count (single-strand families
    that have no duplex partner).
    """
    parsed: Dict[str, List[Dict]] = {}

    for f in module.find_log_files("fgbio/collectduplexseqmetrics_familysizes"):
        module.add_data_source(f)
        s_name = f["s_name"]
        rows: List[Dict] = []
        for line in f["f"].splitlines():
            if line.startswith("family_size"):
                continue
            fields = line.split("\t")
            if len(fields) < 10:
                continue
            family_size = int(fields[0])
            ss_count = int(float(fields[4]))
            ds_count = int(float(fields[7]))
            orphan_count = max(0, ss_count - ds_count)
            rows.append({
                "family_size": family_size,
                "ss_count": ss_count,
                "ds_count": ds_count,
                "orphan_count": orphan_count,
            })

        if rows:
            parsed[s_name] = rows
        module.add_software_version(None, s_name)

    parsed = module.ignore_samples(parsed)
    return parsed


def _add_family_size_plot(module: BaseMultiqcModule, data: Dict[str, List[Dict]]) -> None:
    """Add a line graph showing raw read fractions by family size.

    Three tabs: duplex fraction, orphaned fraction, and raw SSCS family counts.
    """
    duplex_frac_data = {}
    orphan_frac_data = {}
    duplex_count_data = {}
    orphan_count_data = {}

    for s_name, rows in data.items():
        total_reads = sum(
            (r["ds_count"] + r["orphan_count"]) * r["family_size"]
            for r in rows
        )
        if total_reads == 0:
            continue

        d_frac = {}
        o_frac = {}
        d_count = {}
        o_count = {}
        for row in rows:
            fs = row["family_size"]
            if fs > MAX_FAMILY_SIZE:
                continue
            duplex_reads = row["ds_count"] * fs
            orphan_reads = row["orphan_count"] * fs
            d_frac[fs] = duplex_reads / total_reads
            o_frac[fs] = orphan_reads / total_reads
            d_count[fs] = row["ds_count"]
            o_count[fs] = row["orphan_count"]

        duplex_frac_data[f"{s_name} - Duplex"] = d_frac
        orphan_frac_data[f"{s_name} - Orphan"] = o_frac
        duplex_count_data[f"{s_name} - Duplex"] = d_count
        orphan_count_data[f"{s_name} - Orphan"] = o_count

    frac_data = {**duplex_frac_data, **orphan_frac_data}
    count_data = {**duplex_count_data, **orphan_count_data}

    pconfig = {
        "id": "fgbio-duplex-family-sizes-plot",
        "title": "fgbio: SSCS Family Size Distribution",
        "xlab": "Family Size (raw reads per SSCS)",
        "xmin": 1,
        "xmax": MAX_FAMILY_SIZE,
        "x_decimals": False,
        "colors": {},
        "data_labels": [
            {
                "name": "Read Fractions",
                "ylab": "Fraction of Total Raw Reads",
            },
            {
                "name": "Family Counts",
                "ylab": "Number of SSCS Families",
            },
        ],
    }

    for s_name in data:
        pconfig["colors"][f"{s_name} - Duplex"] = "#1b3a5c"
        pconfig["colors"][f"{s_name} - Orphan"] = "#5ec4b5"

    module.add_section(
        name="CollectDuplexSeqMetrics: Family Sizes",
        anchor="fgbio-collectduplexseqmetrics-familysizes",
        description=(
            "Distribution of raw reads across SSCS family sizes from "
            "<code>CollectDuplexSeqMetrics</code>. Duplex SSCS are part of a "
            "complete duplex consensus; orphaned SSCS have no duplex partner."
        ),
        helptext="""
        This plot shows the fraction of total raw reads contributed by each
        single-strand consensus sequence (SSCS) family size. Each sample
        produces two lines:

        - **Duplex** (dark navy): families where both strands (AB and BA)
          are present, forming a complete duplex consensus.
        - **Orphaned** (teal): families where only one strand is present,
          so no duplex consensus can be formed.

        The family size is the number of raw reads that were grouped together
        to form a single SSCS. Higher family sizes indicate deeper per-strand
        coverage. A large fraction of orphaned SSCS at low family sizes (1-2)
        is typical and expected.
        """,
        plot=linegraph.plot([frac_data, count_data], pconfig),
    )

    module.write_data_file(
        {s: {r["family_size"]: r for r in rows} for s, rows in data.items()},
        "fgbio_duplex_family_sizes",
    )


def _add_general_stats(module: BaseMultiqcModule, data: Dict[str, List[Dict]]) -> None:
    """Add duplex rate and mean family size to the general stats table."""
    stats = {}
    for s_name, rows in data.items():
        total_ss = sum(r["ss_count"] for r in rows)
        total_ds = sum(r["ds_count"] for r in rows)
        duplex_rate = total_ds / total_ss if total_ss > 0 else 0.0

        weighted_sum = sum(r["ss_count"] * r["family_size"] for r in rows)
        mean_family_size = weighted_sum / total_ss if total_ss > 0 else 0.0

        stats[s_name] = {
            "duplex_rate": duplex_rate,
            "mean_family_size": mean_family_size,
        }

    headers = {
        "duplex_rate": {
            "title": "Duplex Rate",
            "description": "Fraction of SSCS families that are part of a duplex consensus",
            "min": 0,
            "max": 1,
            "scale": "RdYlGn",
            "format": "{:,.3f}",
        },
        "mean_family_size": {
            "title": "Mean Fam. Size",
            "description": "Mean SSCS family size (weighted by family count)",
            "min": 0,
            "scale": "Blues",
            "format": "{:,.1f}",
        },
    }
    module.general_stats_addcols(stats, headers)
