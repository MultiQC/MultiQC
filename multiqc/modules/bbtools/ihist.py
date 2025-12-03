"""MultiQC submodule to parse BBTools ihist output (insert size histogram)"""

import logging
from itertools import chain
from typing import Dict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, table

log = logging.getLogger(__name__)


def parse_bbtools_ihist(module: BaseMultiqcModule) -> int:
    """
    Parse BBTools ihist output files showing insert size histogram.

    Shows the distribution of computed insert sizes for paired reads.
    """

    data_by_sample: Dict[str, Dict] = {}

    for f in module.find_log_files("bbtools/ihist", filehandles=True):
        parsed = _parse_ihist_file(f["f"])
        if parsed:
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="ihist")
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
        "id": "bbtools-ihist-plot",
        "title": "BBTools: Insert Sizes",
        "xlab": "Insert size (base pairs)",
        "xsuffix": " bp",
        "ylab": "Read count",
    }

    module.add_section(
        name="Insert Sizes",
        anchor="bbtools-ihist",
        description="Histogram of computed insert sizes, for paired reads (`ihist`). "
        "Plotted data has been cut off at 99% to prevent long tails; "
        "complete data available in original source files.",
        helptext="The insert size is the length of the sequence between the "
        "sequencing adapters, which for most common insert sizes is "
        "longer than the sum of both read pairs. "
        "In some cases, the insert size is shorter than the length "
        "of the two read pairs combined, resulting in an insert size "
        "shorter than the sum of the length of the reads pairs.",
        plot=linegraph.plot(plot_data, pconfig),
    )

    # Create summary table if we have key-value data
    table_data = {s: d.get("kv", {}) for s, d in data_by_sample.items() if d.get("kv")}
    if table_data and any(table_data.values()):
        table_headers: Dict = {
            "Mean": {
                "title": "Mean",
                "description": "Average insert length",
            },
            "Median": {
                "title": "Median",
                "description": "Median insert length",
            },
            "STDev": {
                "title": "Std Dev",
                "description": "Standard deviation of insert size length distribution",
            },
            "PercentOfPairs": {
                "title": "% of Pairs",
                "description": "Percentage of pairs",
            },
        }

        module.add_section(
            name="Insert Sizes Summary",
            anchor="bbtools-ihist-table",
            description="Insert size statistics summary.",
            plot=table.plot(
                table_data,
                table_headers,
                pconfig={
                    "id": "bbtools-ihist-summary-table",
                    "namespace": "BBTools",
                    "title": "BBTools: Insert Size Summary",
                },
            ),
        )

    # Write parsed data to file
    module.write_data_file(data_by_sample, "multiqc_bbtools_ihist")

    return len(data_by_sample)


def _parse_ihist_file(f) -> Dict:
    """
    Parse a BBTools ihist file.

    Expected columns: #InsertSize, Count
    Also contains key-value pairs like Mean, Median, STDev, PercentOfPairs
    """
    parsed_data: Dict = {"data": {}, "kv": {}}
    kv_keys = ["Mean", "Median", "STDev", "PercentOfPairs"]

    for line in f:
        line = line.strip()
        if not line:
            continue

        parts = line.split("\t")

        # Check for header/kv lines
        if parts[0].startswith("#"):
            parts[0] = parts[0][1:]  # Remove leading '#'

            # Skip header
            if parts[0] == "InsertSize":
                continue

            # Key-value pair
            if len(parts) == 2 and parts[0] in kv_keys:
                try:
                    parsed_data["kv"][parts[0]] = float(parts[1])
                except ValueError:
                    parsed_data["kv"][parts[0]] = parts[1]
                continue

        # Data row
        if len(parts) >= 2:
            try:
                insert_size = int(parts[0].lstrip("#"))
                count = int(parts[1])
                parsed_data["data"][insert_size] = [count]
            except (ValueError, IndexError):
                continue

    return parsed_data if parsed_data["data"] else {}


def _prepare_plot_data(data_by_sample: Dict) -> Dict:
    """Prepare data for line graph plotting with 99% cutoff to avoid long tails."""
    # Calculate xmax based on 99% of data (not 99.9% like other histograms)
    sumy = sum(
        int(data_by_sample[sample]["data"][x][0])
        for sample in data_by_sample
        for x in data_by_sample[sample].get("data", {})
    )

    cutoff = sumy * 0.99
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
