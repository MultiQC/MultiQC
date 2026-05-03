"""Parse riker `alignment` (alignment-metrics.txt) outputs."""

import logging
from typing import Dict

from multiqc.plots import bargraph, table
from multiqc.plots.bargraph import BarPlotConfig
from multiqc.plots.table import TableConfig

from .util import read_tsv, to_float, to_int

log = logging.getLogger(__name__)

# Column names that should be parsed as integers; everything else is parsed as float
# (read in as strings first so we tolerate the occasional empty/NaN cell).
_INT_COLS = {
    "total_reads",
    "aligned_reads",
    "hq_aligned_reads",
    "hq_aligned_bases",
    "hq_aligned_q20_bases",
    "min_read_length",
    "max_read_length",
    "aligned_reads_in_pairs",
    "reads_improperly_paired",
    "bad_cycles",
}


def parse_reports(module):
    # data_by_sample[sample][category] -> dict of metrics
    data_by_sample: Dict[str, Dict[str, Dict[str, float]]] = {}

    for f in module.find_log_files("riker/alignment", filehandles=True):
        for row in read_tsv(f["f"]):
            sample = row.pop("sample", None)
            category = row.pop("category", None)
            if not sample or not category:
                continue
            s_name = module.clean_s_name(sample, f)

            parsed: Dict[str, float] = {}
            for col, val in row.items():
                if col in _INT_COLS:
                    try:
                        parsed[col] = to_int(val)
                    except (TypeError, ValueError):
                        parsed[col] = to_float(val)
                else:
                    parsed[col] = to_float(val)

            data_by_sample.setdefault(s_name, {})[category] = parsed
            module.add_data_source(f, s_name, section="alignment")

    data_by_sample = module.ignore_samples(data_by_sample)
    if not data_by_sample:
        return set()

    module.add_software_version(None)

    # Pull the `pair` row (or fall back to the first available category) for general stats.
    pair_data: Dict[str, Dict[str, float]] = {}
    for s_name, by_cat in data_by_sample.items():
        pair_data[s_name] = by_cat.get("pair", next(iter(by_cat.values())))

    # Coloring rationale:
    #   - frac_aligned: scale tightly to the meaningful range (80-100%) so
    #     small differences between e.g. 99.5 and 99.9 are visible.
    #   - mismatch_rate / indel_rate / frac_chimeras: lower is better.
    #     Cap the gradient at a "concerning" value so typical-good samples
    #     don't all paint the same shade.
    #   - mean_read_length / total_reads: experiment-dependent. Single-hue
    #     GnBu so the cell shows magnitude without implying good/bad.
    headers = {
        "total_reads": {
            "title": "Total reads",
            "description": "Total number of reads (including QC-failed)",
            "min": 0,
            "format": "{:,.0f}",
            "scale": "GnBu",
            "hidden": True,
        },
        "frac_aligned": {
            "title": "% Aligned",
            "description": "Fraction of PF reads that aligned to the reference",
            "min": 80,
            "max": 100,
            "suffix": "%",
            "format": "{:,.1f}",
            "scale": "RdYlGn",
            "modify": lambda x: x * 100.0,
        },
        "mismatch_rate": {
            "title": "Mismatch rate",
            "description": "Mismatches per aligned base across all PF aligned reads",
            "min": 0,
            "max": 0.015,
            "format": "{:,.4f}",
            "scale": "OrRd",
            "hidden": True,
        },
        "indel_rate": {
            "title": "Indel rate",
            "description": "Insertion + deletion events per aligned base",
            "min": 0,
            "max": 0.005,
            "format": "{:,.4f}",
            "scale": "OrRd",
            "hidden": True,
        },
        "frac_chimeras": {
            "title": "% Chimeric",
            "description": "Fraction of read pairs that are chimeric",
            "min": 0,
            "max": 10,
            "suffix": "%",
            "format": "{:,.2f}",
            "scale": "OrRd",
            "modify": lambda x: x * 100.0,
            "hidden": True,
        },
        "mean_read_length": {
            "title": "Read length",
            "description": "Mean read length across all PF reads",
            "min": 0,
            "suffix": " bp",
            "format": "{:,.0f}",
            "scale": "GnBu",
            "hidden": True,
        },
    }
    module.general_stats_addcols(pair_data, headers, namespace="alignment")

    # Bar plot: aligned vs unaligned reads using the `pair` row.
    bar_data: Dict[str, Dict[str, float]] = {}
    for s_name, by_cat in data_by_sample.items():
        row = by_cat.get("pair", next(iter(by_cat.values())))
        total = row.get("total_reads", 0.0)
        aligned = row.get("aligned_reads", 0.0)
        bar_data[s_name] = {
            "aligned_reads": aligned,
            "unaligned_reads": max(total - aligned, 0.0),
        }

    bar_keys = {
        "aligned_reads": {"name": "Aligned reads"},
        "unaligned_reads": {"name": "Unaligned reads"},
    }
    bar_config = BarPlotConfig(
        id=f"{module.anchor}_alignment_summary",
        title=f"{module.name}: Alignment summary",
        ylab="# Reads",
    )
    module.add_section(
        name="Alignment summary",
        anchor=f"{module.anchor}-alignment-summary",
        description="Aligned vs unaligned reads (from the `pair` row of riker's alignment metrics).",
        plot=bargraph.plot(bar_data, bar_keys, bar_config),
    )

    # Per-category metrics table
    table_data: Dict[str, Dict[str, float]] = {}
    for s_name, by_cat in data_by_sample.items():
        for category, row in by_cat.items():
            table_data[f"{s_name} ({category})"] = row

    table_headers = {
        "total_reads": {"title": "Total reads", "format": "{:,.0f}", "min": 0},
        "aligned_reads": {"title": "Aligned reads", "format": "{:,.0f}", "min": 0},
        "frac_aligned": {
            "title": "% Aligned",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "format": "{:,.2f}",
            "modify": lambda x: x * 100.0,
        },
        "mean_read_length": {"title": "Mean read length", "format": "{:,.1f}", "suffix": " bp"},
        "mismatch_rate": {"title": "Mismatch rate", "format": "{:,.4f}"},
        "indel_rate": {"title": "Indel rate", "format": "{:,.4f}"},
        "frac_chimeras": {
            "title": "% Chimeric",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "format": "{:,.2f}",
            "modify": lambda x: x * 100.0,
        },
        "strand_balance": {"title": "Strand balance", "format": "{:,.3f}"},
        "bad_cycles": {"title": "Bad cycles", "format": "{:,.0f}"},
    }
    table_config = TableConfig(
        id=f"{module.anchor}_alignment_table",
        title=f"{module.name}: Alignment metrics by read category",
    )
    module.add_section(
        name="Alignment metrics",
        anchor=f"{module.anchor}-alignment-table",
        description="Per-category alignment metrics from riker's `alignment` tool.",
        plot=table.plot(table_data, table_headers, table_config),
    )

    module.write_data_file(pair_data, f"multiqc_{module.anchor}_alignment")

    return set(data_by_sample.keys())
