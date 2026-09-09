"""Parse the clipping metrics that fgbio ClipBam writes with `--metrics`."""

import logging
from typing import Dict, List, Union

from multiqc import config
from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table
from multiqc.plots.bargraph import BarPlotConfig
from multiqc.plots.table import TableConfig
from multiqc.plots.table_object import ColumnDict, InputRow, ValueT
from multiqc.types import ColumnKey, SampleGroup, SampleName

log = logging.getLogger(__name__)

# Every column after `read_type` is an integer count.
COUNT_COLUMNS = (
    "reads",
    "reads_unmapped",
    "reads_clipped_pre",
    "reads_clipped_post",
    "reads_clipped_five_prime",
    "reads_clipped_three_prime",
    "reads_clipped_overlapping",
    "reads_clipped_extending",
    "bases",
    "bases_clipped_pre",
    "bases_clipped_post",
    "bases_clipped_five_prime",
    "bases_clipped_three_prime",
    "bases_clipped_overlapping",
    "bases_clipped_extending",
)

# `Pair` sums `ReadOne` and `ReadTwo`, and `All` sums `Fragment` and `Pair`, so `All`
# heads each sample's table group and the other read types nest under it in this order.
HEADLINE_READ_TYPE = "All"
NESTED_READ_TYPES = ("Pair", "ReadOne", "ReadTwo", "Fragment")

# The clipping reasons whose base counts add up to `bases_clipped_post`.
CLIPPED_BASES_BY_REASON = (
    "bases_clipped_pre",
    "bases_clipped_five_prime",
    "bases_clipped_three_prime",
    "bases_clipped_overlapping",
    "bases_clipped_extending",
)

BASES_DENOMINATOR_NOTE = (
    "`bases` in the metrics file counts the aligned bases that remain after clipping, so the "
    "bases present before ClipBam ran are `bases + bases_clipped_post`, which is the denominator "
    "of every base percentage here."
)

ClipBamRow = Dict[str, Union[int, float]]


def run_clip_bam(module: BaseMultiqcModule) -> int:
    """
    Parse fgbio ClipBam metrics files, then add general stats columns, a bar plot of
    clipped bases by reason, and a per-read-type table grouped by sample.
    """
    data_by_sample = parse_clip_bam_metrics(module)
    if not data_by_sample:
        return 0

    module.write_data_file(flatten_by_read_type(data_by_sample), "multiqc_fgbio_clipbam")
    add_general_stats(module, data_by_sample)
    add_clipped_bases_plot(module, data_by_sample)
    add_read_type_table(module, data_by_sample)

    return len(data_by_sample)


def parse_clip_bam_metrics(module: BaseMultiqcModule) -> Dict[str, Dict[str, ClipBamRow]]:
    """Return `{sample: {read_type: metrics}}` for every ClipBam metrics file found."""
    data_by_sample: Dict[str, Dict[str, ClipBamRow]] = {}

    for f in module.find_log_files("fgbio/clipbam"):
        lines = [line for line in f["f"].splitlines() if line.strip()]
        if not lines:
            continue
        header = lines[0].split("\t")
        if header[0] != "read_type" or any(col not in header for col in COUNT_COLUMNS):
            log.debug(f"Skipping {f['fn']}: not a ClipBam metrics header")
            continue

        rows_by_read_type: Dict[str, ClipBamRow] = {}
        for line_num, line in enumerate(lines[1:], start=2):
            fields = line.split("\t")
            if len(fields) != len(header):
                log.warning(
                    f"fgbio ClipBam: skipping line {line_num} of {f['fn']}: "
                    f"got {len(fields)} fields, expected {len(header)}"
                )
                continue
            record = dict(zip(header, fields))
            try:
                counts = {col: int(record[col]) for col in COUNT_COLUMNS}
            except ValueError:
                log.warning(f"fgbio ClipBam: skipping line {line_num} of {f['fn']}: non-integer count")
                continue
            rows_by_read_type[record["read_type"]] = with_derived_metrics(counts)

        if not rows_by_read_type:
            continue
        if f["s_name"] in data_by_sample:
            log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
        data_by_sample[f["s_name"]] = rows_by_read_type
        module.add_data_source(f, section="ClipBam")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        module.add_software_version(None, f["s_name"])

    return module.ignore_samples(data_by_sample)


def with_derived_metrics(counts: Dict[str, int]) -> ClipBamRow:
    """
    Add the percentages the report shows to a row of raw counts.

    `bases` is the number of aligned bases left after clipping, not the read length, so the
    bases present before ClipBam ran are `bases + bases_clipped_post` and every base
    percentage divides by that sum. Read percentages divide by `reads`. A percentage whose
    denominator is zero is left out of the row instead of being reported as zero.
    """
    row: ClipBamRow = dict(counts)
    bases_pre_clip = counts["bases"] + counts["bases_clipped_post"]
    row["bases_pre_clip"] = bases_pre_clip
    if counts["reads"] > 0:
        row["pct_reads_clipped"] = 100.0 * counts["reads_clipped_post"] / counts["reads"]
    if bases_pre_clip > 0:
        row["pct_bases_clipped"] = 100.0 * counts["bases_clipped_post"] / bases_pre_clip
        row["pct_bases_clipped_overlapping"] = 100.0 * counts["bases_clipped_overlapping"] / bases_pre_clip
    return row


def flatten_by_read_type(data_by_sample: Dict[str, Dict[str, ClipBamRow]]) -> Dict[str, ClipBamRow]:
    """One row per sample for the data file, with columns named `<read_type>_<metric>`."""
    return {
        s_name: {
            f"{read_type}_{metric}": value for read_type, row in by_read_type.items() for metric, value in row.items()
        }
        for s_name, by_read_type in data_by_sample.items()
    }


def add_general_stats(module: BaseMultiqcModule, data_by_sample: Dict[str, Dict[str, ClipBamRow]]) -> None:
    gs_keys = ("pct_bases_clipped", "pct_reads_clipped")
    gs_data: Dict[Union[SampleName, str], Dict[Union[ColumnKey, str], ValueT]] = {}
    for s_name, by_read_type in data_by_sample.items():
        all_row = by_read_type.get(HEADLINE_READ_TYPE, {})
        values: Dict[Union[ColumnKey, str], ValueT] = {key: all_row[key] for key in gs_keys if key in all_row}
        if values:
            gs_data[s_name] = values
    if not gs_data:
        return

    headers: Dict[str, ColumnDict] = {
        "pct_bases_clipped": {
            "title": "% Bases clipped",
            "description": (
                "Percentage of bases clipped after ClipBam ran, including clipping already present in the "
                "input: bases_clipped_post / (bases + bases_clipped_post), where bases counts only the aligned "
                "bases left after clipping (All reads)"
            ),
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "OrRd",
            "format": "{:,.1f}",
        },
        "pct_reads_clipped": {
            "title": "% Reads clipped",
            "description": (
                "Percentage of reads with any clipping after ClipBam ran, including clipping already present "
                "in the input: reads_clipped_post / reads (All reads)"
            ),
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "YlOrBr",
            "format": "{:,.1f}",
            "hidden": True,
        },
    }
    module.general_stats_addcols(gs_data, headers, namespace="ClipBam")


def add_clipped_bases_plot(module: BaseMultiqcModule, data_by_sample: Dict[str, Dict[str, ClipBamRow]]) -> None:
    categories = {
        "bases": {"name": "Aligned, not clipped"},
        "bases_clipped_pre": {"name": "Clipped before ClipBam"},
        "bases_clipped_five_prime": {"name": "Clipped at 5' end"},
        "bases_clipped_three_prime": {"name": "Clipped at 3' end"},
        "bases_clipped_overlapping": {"name": "Clipped for mate overlap"},
        "bases_clipped_extending": {"name": "Clipped past mate end"},
    }
    plot_data: Dict[str, Dict[str, Union[int, float]]] = {}
    for s_name, by_read_type in data_by_sample.items():
        all_row = by_read_type.get(HEADLINE_READ_TYPE)
        if all_row is None or all_row["bases_pre_clip"] == 0:
            continue
        plot_data[s_name] = {key: all_row[key] for key in categories}
    if not plot_data:
        return

    module.add_section(
        name="ClipBam clipped bases",
        anchor="fgbio-clipbam-bases",
        description=(
            "Bases in each sample's `All` reads, split into the aligned bases ClipBam left in place and the "
            "bases clipped for each reason."
        ),
        helptext=f"""
        Each bar totals the bases present before ClipBam ran ({BASES_DENOMINATOR_NOTE}).
        The clipped segments are the reasons ClipBam records for a clipped base: clipping that was
        already in the input, the fixed number of bases requested from the 5' or 3' end of a read,
        the bases of a read that overlap its mate, and the bases that extend past the far end of the
        mate. These five counts add up to `bases_clipped_post`.
        """,
        plot=bargraph.plot(
            plot_data,
            categories,
            BarPlotConfig(
                id="fgbio-clipbam-bases-plot",
                title="fgbio: ClipBam bases by clipping reason",
                ylab="# Bases",
                cpswitch_counts_label="Number of bases",
                cpswitch_percent_label="Percentage of bases",
            ),
        ),
    )


def add_read_type_table(module: BaseMultiqcModule, data_by_sample: Dict[str, Dict[str, ClipBamRow]]) -> None:
    # The first row of each group is the headline: `All`, labelled with the bare sample name.
    # The other read types nest under it and expand on click. Read types with no reads are
    # dropped, and a single remaining nested read type is dropped too, since fragment-only
    # data has `Fragment` equal to `All` and the row would only repeat the headline.
    rows_by_sample: Dict[Union[str, SampleGroup], List[InputRow]] = {}
    for s_name, by_read_type in data_by_sample.items():
        headline = by_read_type.get(HEADLINE_READ_TYPE)
        if headline is None or headline["reads"] == 0:
            continue
        nested = [rt for rt in NESTED_READ_TYPES if by_read_type.get(rt, {}).get("reads", 0) > 0]
        if len(nested) == 1:
            nested = []
        rows = [InputRow(sample=SampleName(s_name), data=headline)]
        rows += [InputRow(sample=SampleName(f"{s_name} ({rt})"), data=by_read_type[rt]) for rt in nested]
        rows_by_sample[SampleGroup(s_name)] = rows
    if not rows_by_sample:
        return

    headers: Dict[Union[str, ColumnKey], ColumnDict] = {
        "reads": {
            "title": f"{config.read_count_prefix} Reads",
            "description": f"Reads examined by ClipBam ({config.read_count_desc})",
            "min": 0,
            "format": "{:,.2f}",
            "scale": "GnBu",
            "shared_key": "read_count",
            "modify": lambda x: float(x) * config.read_count_multiplier,
            "hidden": True,
        },
        "pct_reads_clipped": {
            "title": "% Reads clipped",
            "description": (
                "Reads with any clipping after ClipBam ran, including clipping already present in the input: "
                "reads_clipped_post / reads"
            ),
            "min": 0,
            "max": 100,
            "suffix": "%",
            "format": "{:,.1f}",
            "scale": "YlOrBr",
        },
        "bases_pre_clip": {
            "title": f"{config.base_count_prefix} Bases pre-clip",
            "description": (
                f"Bases present before ClipBam ran: bases + bases_clipped_post, where bases counts only the "
                f"aligned bases left after clipping ({config.base_count_desc})"
            ),
            "min": 0,
            "format": "{:,.2f}",
            "scale": "Blues",
            "shared_key": "base_count",
            "modify": lambda x: float(x) * config.base_count_multiplier,
            "hidden": True,
        },
        "pct_bases_clipped": {
            "title": "% Bases clipped",
            "description": (
                "Bases clipped after ClipBam ran, including clipping already present in the input: "
                "bases_clipped_post / (bases + bases_clipped_post)"
            ),
            "min": 0,
            "max": 100,
            "suffix": "%",
            "format": "{:,.1f}",
            "scale": "OrRd",
        },
        "pct_bases_clipped_overlapping": {
            "title": "% Bases clipped (overlap)",
            "description": (
                "Bases clipped because the read overlapped its mate: "
                "bases_clipped_overlapping / (bases + bases_clipped_post)"
            ),
            "min": 0,
            "max": 100,
            "suffix": "%",
            "format": "{:,.1f}",
            "scale": "PuRd",
        },
        "bases_clipped_overlapping": {
            "title": f"{config.base_count_prefix} Clipped (overlap)",
            "description": f"Bases clipped because the read overlapped its mate ({config.base_count_desc})",
            "min": 0,
            "format": "{:,.2f}",
            "scale": "Purples",
            "shared_key": "base_count",
            "modify": lambda x: float(x) * config.base_count_multiplier,
        },
        "bases_clipped_pre": {
            "title": f"{config.base_count_prefix} Clipped (pre-existing)",
            "description": f"Bases already clipped in the input before ClipBam ran ({config.base_count_desc})",
            "min": 0,
            "format": "{:,.2f}",
            "scale": "Greys",
            "shared_key": "base_count",
            "modify": lambda x: float(x) * config.base_count_multiplier,
            "hidden": True,
        },
        "bases_clipped_five_prime": {
            "title": f"{config.base_count_prefix} Clipped (5')",
            "description": f"Bases clipped from the 5' end of reads ({config.base_count_desc})",
            "min": 0,
            "format": "{:,.2f}",
            "scale": "Greens",
            "shared_key": "base_count",
            "modify": lambda x: float(x) * config.base_count_multiplier,
            "hidden": True,
        },
        "bases_clipped_three_prime": {
            "title": f"{config.base_count_prefix} Clipped (3')",
            "description": f"Bases clipped from the 3' end of reads ({config.base_count_desc})",
            "min": 0,
            "format": "{:,.2f}",
            "scale": "YlGn",
            "shared_key": "base_count",
            "modify": lambda x: float(x) * config.base_count_multiplier,
            "hidden": True,
        },
        "bases_clipped_extending": {
            "title": f"{config.base_count_prefix} Clipped (past mate)",
            "description": f"Bases clipped because the read extended past the far end of its mate ({config.base_count_desc})",
            "min": 0,
            "format": "{:,.2f}",
            "scale": "Oranges",
            "shared_key": "base_count",
            "modify": lambda x: float(x) * config.base_count_multiplier,
            "hidden": True,
        },
    }

    module.add_section(
        name="ClipBam clipping metrics",
        anchor="fgbio-clipbam",
        description=(
            "Reads and bases clipped by `ClipBam`, per sample and read type. Each sample's headline row "
            "holds the `All` metrics; expand it for the `Pair`, `ReadOne`, `ReadTwo` and `Fragment` rows."
        ),
        helptext=f"""
        `ClipBam` clips reads in an aligned BAM: a fixed number of bases from the 5' or 3' end,
        the bases of a read that overlap its mate, or the bases that extend past the far end of
        the mate. Its metrics file counts, for each read type, how many reads and bases ended up
        clipped and for which reason.

        Read types: `ReadOne` and `ReadTwo` are the two reads of a pair, `Pair` is their sum,
        `Fragment` is unpaired reads, and `All` is `Pair` plus `Fragment`. The `All` row is each
        sample's headline row; the other read types nest under it and can be shown by clicking
        the row. Read types with no reads are not shown, so paired-end data has no `Fragment` row
        and fragment-only data collapses to the headline row alone.

        Denominators: {BASES_DENOMINATOR_NOTE}
        `% Bases clipped` is `bases_clipped_post / (bases + bases_clipped_post)` and
        `% Bases clipped (overlap)` is `bases_clipped_overlapping / (bases + bases_clipped_post)`.
        `% Reads clipped` is `reads_clipped_post / reads`. The `_post` counts include clipping
        that was already present in the input (the `_pre` counts), and the per-reason base
        counts (pre-existing, 5', 3', overlap, past mate) add up to `bases_clipped_post`.
        """,
        plot=table.plot(
            rows_by_sample,
            headers,
            TableConfig(
                id="fgbio-clipbam-table",
                title="fgbio: ClipBam clipping metrics by read type",
            ),
        ),
    )
