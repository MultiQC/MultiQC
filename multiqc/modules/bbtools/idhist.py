"""MultiQC submodule to parse BBTools idhist output (identity histogram)"""

import logging
from itertools import chain
from typing import Dict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, table

log = logging.getLogger(__name__)


def parse_bbtools_idhist(module: BaseMultiqcModule) -> int:
    """
    Parse BBTools idhist output files showing histogram of read identity percentages.

    Shows read count versus percent base pair identity of aligned reads.
    """

    data_by_sample: Dict[str, Dict] = {}

    for f in module.find_log_files("bbtools/idhist", filehandles=True):
        parsed = _parse_idhist_file(f["f"])
        if parsed:
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="idhist")
            data_by_sample[s_name] = parsed

    # Filter to strip out ignored sample names
    data_by_sample = module.ignore_samples(data_by_sample)

    if len(data_by_sample) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    module.add_software_version(None)

    # Create line graph plot
    plot_data = _prepare_plot_data(data_by_sample)

    pconfig = {
        "id": "bbtools-idhist-plot",
        "title": "BBTools: Identity Histogram",
        "xlab": "Percent identity",
        "xsuffix": "%",
        "ylab": "Read count",
        "data_labels": [
            {"name": "Reads", "ylab": "Read count"},
            {"name": "Bases", "ylab": "Number of bases"},
        ],
    }

    module.add_section(
        name="Identity Histogram",
        anchor="bbtools-idhist",
        description="Histogram of read count versus percent base pair identity of aligned reads (`idhist`).",
        plot=linegraph.plot(plot_data, pconfig),
    )

    # Create summary table if we have key-value data
    table_data = {s: d.get("kv", {}) for s, d in data_by_sample.items() if d.get("kv")}
    if table_data and any(table_data.values()):
        table_headers: Dict = {
            "Mean_reads": {
                "title": "Mean (reads)",
                "description": "Average percent identity of aligned reads",
            },
            "Mean_bases": {
                "title": "Mean (bases)",
                "description": "Average percent identity of aligned bases",
            },
            "Median_reads": {
                "title": "Median (reads)",
                "description": "Median percent identity of aligned reads",
            },
            "Median_bases": {
                "title": "Median (bases)",
                "description": "Median percent identity of aligned bases",
            },
            "Mode_reads": {
                "title": "Mode (reads)",
                "description": "The most commonly occurring value amongst aligned reads",
            },
            "Mode_bases": {
                "title": "Mode (bases)",
                "description": "The most commonly occurring value amongst aligned bases",
            },
            "STDev_reads": {
                "title": "Std Dev (reads)",
                "description": "Standard deviation for aligned reads",
            },
            "STDev_bases": {
                "title": "Std Dev (bases)",
                "description": "Standard deviation for aligned bases",
            },
        }

        module.add_section(
            name="Identity Histogram Summary",
            anchor="bbtools-idhist-table",
            description="Identity statistics summary.",
            plot=table.plot(
                table_data,
                table_headers,
                pconfig={
                    "id": "bbtools-idhist-summary-table",
                    "namespace": "BBTools",
                    "title": "BBTools: Identity Summary",
                },
            ),
        )

    # Write parsed data to file
    module.write_data_file(data_by_sample, "multiqc_bbtools_idhist")

    return len(data_by_sample)


def _parse_idhist_file(f) -> Dict:
    """
    Parse a BBTools idhist file.

    Expected columns: #Identity, Reads, Bases
    Also contains key-value pairs like Mean_reads, Mean_bases, etc.
    """
    parsed_data: Dict = {"data": {}, "kv": {}}
    kv_keys = [
        "Mean_reads", "Mean_bases", "Median_reads", "Median_bases",
        "Mode_reads", "Mode_bases", "STDev_reads", "STDev_bases",
    ]

    for line in f:
        line = line.strip()
        if not line:
            continue

        parts = line.split("\t")

        # Check for header/kv lines
        if parts[0].startswith("#"):
            parts[0] = parts[0][1:]  # Remove leading '#'

            # Skip header
            if parts[0] == "Identity":
                continue

            # Key-value pair
            if len(parts) == 2 and parts[0] in kv_keys:
                try:
                    parsed_data["kv"][parts[0]] = float(parts[1])
                except ValueError:
                    parsed_data["kv"][parts[0]] = parts[1]
                continue

        # Data row
        if len(parts) >= 3:
            try:
                identity = float(parts[0].lstrip("#"))
                reads = int(parts[1])
                bases = int(parts[2])
                parsed_data["data"][identity] = [reads, bases]
            except (ValueError, IndexError):
                continue

    return parsed_data if parsed_data["data"] else {}


def _prepare_plot_data(data_by_sample: Dict) -> list:
    """Prepare data for line graph plotting with Reads and Bases datasets."""
    # Get all x values
    all_x = set()
    for sample_data in data_by_sample.values():
        all_x.update(sample_data.get("data", {}).keys())

    # Columns: Reads = 0, Bases = 1
    columns_to_plot = {
        "Reads": {0: "Count"},
        "Bases": {1: "Count"},
    }

    plot_data = []
    for column_type in columns_to_plot:
        dataset = {}
        for column, column_name in columns_to_plot[column_type].items():
            for sample in data_by_sample:
                sample_key = f"{sample}.{column_name}"
                dataset[sample_key] = {}
                sample_data = data_by_sample[sample].get("data", {})
                for x in all_x:
                    if x in sample_data:
                        dataset[sample_key][x] = sample_data[x][column]
                    else:
                        dataset[sample_key][x] = 0
        plot_data.append(dataset)

    return plot_data
