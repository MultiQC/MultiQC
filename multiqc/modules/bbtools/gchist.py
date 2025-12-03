"""MultiQC submodule to parse BBTools gchist output (GC content histogram)"""

import logging
from itertools import chain
from typing import Dict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, table

log = logging.getLogger(__name__)


def parse_bbtools_gchist(module: BaseMultiqcModule) -> int:
    """
    Parse BBTools gchist output files showing read GC content histogram.
    """

    data_by_sample: Dict[str, Dict] = {}

    for f in module.find_log_files("bbtools/gchist", filehandles=True):
        parsed = _parse_gchist_file(f["f"])
        if parsed:
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="gchist")
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
        "id": "bbtools-gchist-plot",
        "title": "BBTools: GC Content",
        "xlab": "Proportion GC",
        "xsuffix": "%",
        "xmin": 0,
        "xmax": 100,
        "ylab": "# Reads",
    }

    module.add_section(
        name="GC Content",
        anchor="bbtools-gchist",
        description="Read GC content histogram (`gchist`).",
        plot=linegraph.plot(plot_data, pconfig),
    )

    # Create summary table if we have key-value data
    table_data = {s: d.get("kv", {}) for s, d in data_by_sample.items() if d.get("kv")}
    if table_data and any(table_data.values()):
        table_headers: Dict = {
            "Mean": {
                "title": "Mean GC",
                "description": "Average GC content",
            },
            "Median": {
                "title": "Median GC",
                "description": "Median GC content",
            },
            "Mode": {
                "title": "Mode GC",
                "description": "The most commonly occurring value of the GC content distribution",
            },
            "STDev": {
                "title": "Std Dev",
                "description": "Standard deviation of average GC content",
            },
        }

        module.add_section(
            name="GC Content Summary",
            anchor="bbtools-gchist-table",
            description="Read GC content statistics.",
            plot=table.plot(
                table_data,
                table_headers,
                pconfig={
                    "id": "bbtools-gchist-summary-table",
                    "namespace": "BBTools",
                    "title": "BBTools: GC Content Summary",
                },
            ),
        )

    # Write parsed data to file
    module.write_data_file(data_by_sample, "multiqc_bbtools_gchist")

    return len(data_by_sample)


def _parse_gchist_file(f) -> Dict:
    """
    Parse a BBTools gchist file.

    Expected columns: #GC, Count
    Also contains key-value pairs like Mean, Median, Mode, STDev
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
            if parts[0] == "GC":
                continue

            # Key-value pair (Mean, Median, etc.)
            if len(parts) == 2 and parts[0] in ["Mean", "Median", "Mode", "STDev"]:
                try:
                    parsed_data["kv"][parts[0]] = float(parts[1])
                except ValueError:
                    parsed_data["kv"][parts[0]] = parts[1]
                continue

        # Data row
        if len(parts) >= 2:
            try:
                gc = float(parts[0].lstrip("#"))
                count = int(parts[1])
                parsed_data["data"][gc] = [count]
            except (ValueError, IndexError):
                continue

    return parsed_data if parsed_data["data"] else {}


def _prepare_plot_data(data_by_sample: Dict) -> Dict:
    """Prepare data for line graph plotting."""
    # Calculate xmax based on 99.9% of data
    sumy = sum(
        int(data_by_sample[sample]["data"][x][0])
        for sample in data_by_sample
        for x in data_by_sample[sample].get("data", {})
    )

    cutoff = sumy * 0.999
    all_x = set()
    items = []
    for sample in data_by_sample:
        for x, vals in data_by_sample[sample].get("data", {}).items():
            items.append((x, vals))

    for item in sorted(items):
        all_x.add(item[0])
        cutoff -= item[1][0]
        if cutoff < 0:
            break

    # Build plot data
    plot_data = {}
    for sample in data_by_sample:
        plot_data[sample] = {}
        sample_data = data_by_sample[sample].get("data", {})
        for x in all_x:
            if x in sample_data:
                plot_data[sample][x] = sample_data[x][0]
            else:
                plot_data[sample][x] = 0

    return plot_data
