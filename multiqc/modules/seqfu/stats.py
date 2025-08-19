import logging
from typing import Dict, Any

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, SampleGroupingConfig, ModuleNoSamplesFound
from multiqc.plots import bargraph
from multiqc.types import SampleName, ColumnKey

log = logging.getLogger(__name__)


def parse_seqfu_stats(module: BaseMultiqcModule):
    """Find Seqfu stats logs and parse their data"""
    use_filename = False
    if isinstance(config.use_filename_as_sample_name, list):
        # Check for module anchor
        if module.anchor in config.use_filename_as_sample_name:
            use_filename = True
    elif config.use_filename_as_sample_name is True:
        use_filename = True

    seqfu_stats: Dict[SampleName, Dict[str, Any]] = {}
    for f in module.find_log_files("seqfu/stats", filehandles=True):
        for sample_name, data in parse_file(f, use_filename):
            sample_name = SampleName(module.clean_s_name(sample_name, f=f))

            if sample_name in seqfu_stats:
                log.debug(f"Duplicate sample name found! Overwriting: {sample_name}")

            seqfu_stats[sample_name] = data

            module.add_data_source(f, sample_name)

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            module.add_software_version(None, sample_name)

    # Filter to strip out ignored sample names
    seqfu_stats = module.ignore_samples(seqfu_stats)

    if len(seqfu_stats) == 0:
        raise ModuleNoSamplesFound

    add_general_stats_cols(module, seqfu_stats)

    plot_sequence_lengths(module, seqfu_stats)

    plot_sequence_counts(module, seqfu_stats)

    # Write parsed report data to a file
    module.write_data_file(seqfu_stats, "multiqc_seqfu_stats")

    # Return the number of logs that were found
    return len(seqfu_stats)


def parse_file(f, use_filename=False):
    lines = f["f"].readlines()

    # Check if file is empty
    if len(lines) == 1:
        log.debug("Empty file detected")
        return []

    if len(lines) > 2 and use_filename:
        log.warning(f"File with multiple samples incompatible with option `seqfu_stats_config.use_filename` ({f['f']})")
        return []

    # Extract columns from first line
    cols = [c for c in lines[0].strip().split("\t")]

    # Parse sample(s) data
    for line in lines[1:]:
        values = line.strip().split("\t")
        row = dict(zip(cols, values))

        sample_name = row["File"]

        # Using file name if generated from stdin (File='-') or option `use_filename` selected
        if sample_name == "-" or use_filename:
            sample_name = f["s_name"]

        data = {k: float(v) for k, v in row.items() if k != "File"}
        yield sample_name, data


def add_general_stats_cols(module: BaseMultiqcModule, seqfu_stats: Dict[SampleName, Dict[str, Any]]):
    # Add columns to General Stats Table
    general_stats_headers = get_general_stats_headers()

    cols_to_weighted_average = [
        (ColumnKey("Avg"), ColumnKey("#Seq")),
        (ColumnKey("N50"), ColumnKey("#Seq")),
        (ColumnKey("N75"), ColumnKey("#Seq")),
        (ColumnKey("N90"), ColumnKey("#Seq")),
        (ColumnKey("auN"), ColumnKey("#Seq")),
        (ColumnKey("Min"), ColumnKey("#Seq")),
        (ColumnKey("Max"), ColumnKey("#Seq")),
    ]

    no_gc = not any("%GC" in row for row in seqfu_stats.values())
    if not no_gc:
        cols_to_weighted_average.append((ColumnKey("%GC"), ColumnKey("#Seq")))

    general_stats_grouping_config = SampleGroupingConfig(
        cols_to_sum=[ColumnKey(k) for k in ["#Seq", "Total bp"]], cols_to_weighted_average=cols_to_weighted_average
    )

    module.general_stats_addcols(
        data_by_sample={s_name: {ColumnKey(k): v for k, v in data.items()} for s_name, data in seqfu_stats.items()},
        headers=general_stats_headers,
        namespace="stats",
        group_samples_config=general_stats_grouping_config,
    )


def all_same_length(seqfu_stats: Dict[SampleName, Dict[str, Any]]):
    """Check if all sequences are the same length"""
    lengths: set[float] = set()
    for col in ["Min", "Max"]:
        lengths.update(data[col] for data in seqfu_stats.values())
    return len(lengths) == 1


def plot_sequence_lengths(module: BaseMultiqcModule, seqfu_stats: Dict[SampleName, Dict[str, Any]]):
    """
    Plot sequence length statistics as a bar graph with switches for different stats
    """
    if all_same_length(seqfu_stats):
        log.debug("All samples have sequences of a single length")

        # Show a message if all sequences are the same length
        # code inspired by FastQC module
        length = seqfu_stats[next(iter(seqfu_stats))]["Min"]
        module.add_section(
            name="Sequence lengths",
            anchor="seqfu-stats-lengths",
            content=f'<div class="alert alert-info">All samples have sequences of a single length ({length:,.0f} bp)</div>',
        )

        return

    seqfu_lengths_cols = ["Avg", "N50", "N75", "N90", "auN", "Min", "Max"]

    seqfu_lengths_cols_labels = [
        {"name": "Mean", "ylab": "Mean length"},
        {"name": "N50", "ylab": "N50 length"},
        {"name": "N75", "ylab": "N75 length"},
        {"name": "N90", "ylab": "N90 length"},
        {"name": "auN", "ylab": "auN length"},
        {"name": "Min", "ylab": "Min length"},
        {"name": "Max", "ylab": "Max length"},
    ]

    seqfu_lengths_data = []
    for c in seqfu_lengths_cols:
        seqfu_lengths_data.append({str(s): {c: seqfu_stats[s][c]} for s in seqfu_stats.keys()})

    module.add_section(
        name="Sequence lengths",
        anchor="seqfu-stats-lengths",
        id="seqfu-stats-lengths",
        description="Sequence lengths statistics from `seqfu stats`",
        helptext="""
- Mean: average sequence length
- N50: 50% of sequences are longer than this
- N75: 75% of sequences are longer than this
- N90: 90% of sequences are longer than this
- auN: Area under the Nx sequence length curve
- Min: minimum sequence length
- Max: maximum sequence length
""",
        plot=bargraph.plot(
            seqfu_lengths_data,
            pconfig={
                "id": "seqfu-stats-lengths-barplot",
                "title": "Seqfu stats: Sequence length statistics",
                "ymin": 0,
                "cpswitch": False,
                "ysuffix": " bp",
                "tt_decimals": 0,
                "data_labels": seqfu_lengths_cols_labels,
            },
        ),
    )


def plot_sequence_counts(module: BaseMultiqcModule, seqfu_stats: Dict[SampleName, Dict[str, Any]]):
    """
    Plot sequence count statistics as a bar graph with switches for sequences and bases
    """
    seqfu_counts_cols = ["#Seq", "Total bp"]

    seqfu_counts_cols_labels = [
        {"name": "Sequences", "ylab": "Sequences"},
        {"name": "Bases", "ylab": "Bases"},
    ]

    seqfu_counts_data = []
    for c in seqfu_counts_cols:
        seqfu_counts_data.append({str(s): {c: seqfu_stats[s][c]} for s in seqfu_stats.keys()})

    module.add_section(
        name="Sequence counts",
        anchor="seqfu-stats-counts",
        id="seqfu-stats-counts",
        description="Sequence count statistics from `seqfu stats`",
        plot=bargraph.plot(
            seqfu_counts_data,
            pconfig={
                "id": "seqfu-stats-counts-barplot",
                "title": "Seqfu stats: Sequence count statistics",
                "ymin": 0,
                "ysuffix": "",
                "cpswitch": False,
                "tt_decimals": 0,
                "data_labels": seqfu_counts_cols_labels,
            },
        ),
    )


def get_general_stats_headers():
    return {
        "#Seq": {
            "title": "Seqs",
            "description": "Number of sequences",
            "shared_key": "read_count",
            "scale": "Oranges",
        },
        "Total bp": {
            "title": "Bases",
            "description": "Number of bases",
            "shared_key": "base_count",
            "scale": "Purples",
        },
        "Avg": {
            "title": "Mean len",
            "description": "Average sequence length",
            "format": "{:,.0f}",
            "scale": "Greens",
            "suffix": "bp",
        },
        "N50": {
            "title": "N50 len",
            "description": "50% of the sequences are longer than this size",
            "format": "{:,.0f}",
            "scale": "Blues",
            "suffix": "bp",
            "hidden": True,
        },
        "N75": {
            "title": "N75 len",
            "description": "75% of the sequences are longer than this size",
            "format": "{:,.0f}",
            "scale": "Blues",
            "suffix": "bp",
            "hidden": True,
        },
        "N90": {
            "title": "N90 len",
            "description": "90% of the sequences are longer than this size",
            "format": "{:,.0f}",
            "scale": "Blues",
            "suffix": "bp",
            "hidden": True,
        },
        "auN": {
            "title": "auN len",
            "description": "Area under the Nx curve",
            "format": "{:,.0f}",
            "scale": "Blues",
            "suffix": "bp",
            "hidden": True,
        },
        "Min": {
            "title": "Min len",
            "description": "Length of the shortest sequence",
            "format": "{:,.0f}",
            "scale": "RdYlGn",
            "suffix": "bp",
            "hidden": True,
        },
        "Max": {
            "title": "Max len",
            "description": "Length of the longest sequence",
            "format": "{:,.0f}",
            "scale": "RdYlGn",
            "suffix": "bp",
            "hidden": True,
        },
        "%GC": {
            "title": "%GC",
            "description": "GC content in sequences",
            "format": "{:.1%}",
            "scale": "Oranges",
            "min": 0,
            "max": 1,
        },
    }
