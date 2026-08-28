"""Parse riker `alignment` (alignment-metrics.txt) outputs."""

import logging
from typing import Dict

from multiqc import config
from multiqc.plots import bargraph, table
from multiqc.plots.bargraph import BarPlotConfig
from multiqc.plots.table import TableConfig

from .util import read_tsv, to_int

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
        for row in read_tsv(f["f"], source=f["fn"]):
            sample = row.pop("sample", None)
            category = row.pop("category", None)
            if not sample or not category:
                continue
            s_name = module.clean_s_name(sample, f)

            try:
                parsed: Dict[str, float] = {
                    col: (to_int(val) if col in _INT_COLS else float(val)) for col, val in row.items()
                }
            except (TypeError, ValueError) as e:
                log.warning(f"riker: skipping row in {f['fn']} for sample {sample}: {e}")
                continue

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
            "title": f"{config.read_count_prefix} Total reads",
            "description": f"Total number of reads, including QC-failed ({config.read_count_desc})",
            "min": 0,
            "format": "{:,.2f}",
            "scale": "GnBu",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
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
        helptext="""
            Total reads per sample split into aligned and unaligned, taken from the `pair`
            category row of riker's `alignment` output (or the first available category if
            `pair` is not present). The total includes QC-failed reads.
        """,
        plot=bargraph.plot(bar_data, bar_keys, bar_config),
    )

    # Per-category metrics table
    table_data: Dict[str, Dict[str, float]] = {}
    for s_name, by_cat in data_by_sample.items():
        for category, row in by_cat.items():
            table_data[f"{s_name} ({category})"] = row

    # Same coloring rationale as the general stats columns above.
    table_headers = {
        "total_reads": {
            "title": f"{config.read_count_prefix} Total reads",
            "description": f"Total number of reads ({config.read_count_desc})",
            "format": "{:,.2f}",
            "min": 0,
            "scale": "GnBu",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        },
        "aligned_reads": {
            "title": f"{config.read_count_prefix} Aligned reads",
            "description": f"Number of aligned reads ({config.read_count_desc})",
            "format": "{:,.2f}",
            "min": 0,
            "scale": "GnBu",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        },
        "frac_aligned": {
            "title": "% Aligned",
            "min": 80,
            "max": 100,
            "suffix": "%",
            "format": "{:,.2f}",
            "scale": "RdYlGn",
            "modify": lambda x: x * 100.0,
        },
        "mean_read_length": {
            "title": "Mean read length",
            "format": "{:,.1f}",
            "suffix": " bp",
            "min": 0,
            "scale": "GnBu",
        },
        "mismatch_rate": {
            "title": "Mismatch rate",
            "format": "{:,.4f}",
            "min": 0,
            "max": 0.015,
            "scale": "OrRd",
        },
        "indel_rate": {
            "title": "Indel rate",
            "format": "{:,.4f}",
            "min": 0,
            "max": 0.005,
            "scale": "OrRd",
        },
        "frac_chimeras": {
            "title": "% Chimeric",
            "min": 0,
            "max": 10,
            "suffix": "%",
            "format": "{:,.2f}",
            "scale": "OrRd",
            "modify": lambda x: x * 100.0,
        },
        "strand_balance": {
            "title": "Strand balance",
            "format": "{:,.3f}",
            "min": 0,
            "max": 1,
        },
        "bad_cycles": {
            "title": "Bad cycles",
            "format": "{:,.0f}",
            "min": 0,
            "scale": "OrRd",
        },
    }
    table_config = TableConfig(
        id=f"{module.anchor}_alignment_table",
        title=f"{module.name}: Alignment metrics by read category",
    )
    module.add_section(
        name="Alignment metrics",
        anchor=f"{module.anchor}-alignment-table",
        description="Per-category alignment metrics from riker's `alignment` tool.",
        helptext="""
            Alignment summary metrics, equivalent to Picard `CollectAlignmentSummaryMetrics`.
            Each sample is split into rows for the `pair`, `first_of_pair`, and `second_of_pair`
            read categories.

            * `Total reads` / `Aligned reads`: read counts (QC-failed reads are included in the total).
            * `% Aligned`: aligned reads as a fraction of total.
            * `Mean read length`: across all PF (passing filter) reads.
            * `Mismatch rate` / `Indel rate`: per-base mismatch and insertion/deletion rates
              across all aligned bases.
            * `% Chimeric`: pairs that map across contigs, span unexpectedly large inserts,
              or have an unexpected orientation.
            * `Strand balance`: fraction of aligned reads on the forward strand. 0.5 is unbiased.
            * `Bad cycles`: sequencing cycles where at least 80% of reads called N.
        """,
        plot=table.plot(table_data, table_headers, table_config),
    )

    module.write_data_file(pair_data, f"multiqc_{module.anchor}_alignment")

    return set(data_by_sample.keys())
