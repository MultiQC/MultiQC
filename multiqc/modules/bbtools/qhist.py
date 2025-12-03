"""MultiQC submodule to parse BBTools qhist output (quality by position)"""

import logging
from itertools import chain
from typing import Dict, Set

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, table

log = logging.getLogger(__name__)


def parse_bbtools_qhist(module: BaseMultiqcModule) -> int:
    """
    Parse BBTools qhist output files showing quality histogram by position.

    Shows the average quality for each position in the reads.
    """

    data_by_sample: Dict[str, Dict] = {}

    for f in module.find_log_files("bbtools/qhist", filehandles=True):
        parsed = _parse_qhist_file(f["f"])
        if parsed:
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="qhist")
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
        "id": "bbtools-qhist-plot",
        "title": "BBTools: Sequence Quality Histograms",
        "xlab": "Position in read",
        "xsuffix": " bp",
        "ylab": "Quality score",
        "data_labels": [
            {"name": "Linear", "ylab": "Quality score"},
            {"name": "Logarithmic", "ylab": "Log score"},
            {"name": "Measured", "ylab": "Measured quality value"},
        ],
    }

    module.add_section(
        name="Sequence Quality Histograms",
        anchor="bbtools-qhist",
        description="Quality histogram by position (`qhist`). "
        "The plots show the average quality for each position in the reads, "
        "using the linear values, logarithmically scaled values, and the "
        "actual measured qualities based on the alignments.",
        plot=linegraph.plot(plot_data, pconfig),
    )

    # Create summary table if we have key-value data
    table_data = {s: d.get("kv", {}) for s, d in data_by_sample.items() if d.get("kv")}
    if table_data and any(table_data.values()):
        table_headers: Dict = {
            "Deviation": {
                "title": "Deviation",
                "description": "Deviation from expected quality",
            },
        }

        module.add_section(
            name="Quality by Position Summary",
            anchor="bbtools-qhist-table",
            description="Quality by position statistics summary.",
            plot=table.plot(
                table_data,
                table_headers,
                pconfig={
                    "id": "bbtools-qhist-summary-table",
                    "namespace": "BBTools",
                    "title": "BBTools: Quality by Position Summary",
                },
            ),
        )

    # Write parsed data to file
    module.write_data_file(data_by_sample, "multiqc_bbtools_qhist")

    return len(data_by_sample)


def _parse_qhist_file(f) -> Dict:
    """
    Parse a BBTools qhist file.

    Expected columns: #BaseNum, Read1_linear, Read1_log, Read1_measured,
                      [Read2_linear, Read2_log, Read2_measured]
    """
    parsed_data: Dict = {"data": {}, "kv": {}}

    for line in f:
        line = line.strip()
        if not line:
            continue

        parts = line.split("\t")

        # Check for header/kv lines
        if parts[0].startswith("#"):
            parts[0] = parts[0][1:]  # Remove leading '#'

            # Skip header
            if parts[0] == "BaseNum":
                continue

            # Key-value pair (Deviation)
            if len(parts) == 2 and parts[0] == "Deviation":
                try:
                    parsed_data["kv"][parts[0]] = float(parts[1])
                except ValueError:
                    parsed_data["kv"][parts[0]] = parts[1]
                continue

        # Data row
        if len(parts) >= 4:
            try:
                base_num = int(parts[0].lstrip("#"))
                read1_linear = float(parts[1])
                read1_log = float(parts[2])
                read1_measured = float(parts[3])
                # Optional Read2 columns
                read2_linear = float(parts[4]) if len(parts) > 4 else 0.0
                read2_log = float(parts[5]) if len(parts) > 5 else 0.0
                read2_measured = float(parts[6]) if len(parts) > 6 else 0.0
                parsed_data["data"][base_num] = [
                    read1_linear, read1_log, read1_measured,
                    read2_linear, read2_log, read2_measured,
                ]
            except (ValueError, IndexError):
                continue

    return parsed_data if parsed_data["data"] else {}


def _prepare_plot_data(data_by_sample: Dict) -> list:
    """Prepare data for line graph plotting with Linear, Log, and Measured datasets."""
    # Calculate xmax based on 99.9% of data
    sumy = sum(
        data_by_sample[sample]["data"][x][0]
        for sample in data_by_sample
        for x in data_by_sample[sample].get("data", {})
    )

    cutoff = sumy * 0.999
    all_x: Set[int] = set()
    items = []
    for sample in data_by_sample:
        for x, vals in data_by_sample[sample].get("data", {}).items():
            items.append((x, vals))

    for item in sorted(items):
        all_x.add(item[0])
        cutoff -= item[1][0]
        if cutoff < 0:
            break

    # Linear=0, Log=1, Measured=2 for Read1
    # Linear=3, Log=4, Measured=5 for Read2
    columns_to_plot = {
        "Linear": {0: "Read1", 3: "Read2"},
        "Logarithmic": {1: "Read1", 4: "Read2"},
        "Measured": {2: "Read1", 5: "Read2"},
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
                    if x in sample_data and len(sample_data[x]) > column:
                        dataset[sample_key][x] = sample_data[x][column]
        plot_data.append(dataset)

    return plot_data
