"""MultiQC submodule to parse BBTools aqhist output (average read quality histogram)"""

import logging
from itertools import chain
from typing import Dict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


def parse_bbtools_aqhist(module: BaseMultiqcModule) -> int:
    """
    Parse BBTools aqhist output files showing histogram of average read qualities.

    Shows the number of reads at each average quality score.
    """

    data_by_sample: Dict[str, Dict] = {}

    for f in module.find_log_files("bbtools/aqhist", filehandles=True):
        parsed = _parse_aqhist_file(f["f"])
        if parsed:
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="aqhist")
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
        "id": "bbtools-aqhist-plot",
        "title": "BBTools: Average Read Quality",
        "xlab": "Quality score",
        "ylab": "Read count",
        "ylog": True,
        "x_bands": [
            {"from": 28, "to": 100, "color": "#009500", "opacity": 0.13},
            {"from": 20, "to": 28, "color": "#a07300", "opacity": 0.13},
            {"from": 0, "to": 20, "color": "#990101", "opacity": 0.13},
        ],
        "data_labels": [
            {"name": "Count data", "ylab": "Read count"},
            {"name": "Proportion data", "ylab": "Proportion of reads"},
        ],
    }

    module.add_section(
        name="Read Quality",
        anchor="bbtools-aqhist",
        description="Histogram of average read qualities (`aqhist`). "
        "Plot shows the number of reads at each quality score.",
        plot=linegraph.plot(plot_data, pconfig),
    )

    # Write parsed data to file
    module.write_data_file(data_by_sample, "multiqc_bbtools_aqhist")

    return len(data_by_sample)


def _parse_aqhist_file(f) -> Dict:
    """
    Parse a BBTools aqhist file.

    Expected columns: #Quality, count1, fraction1, count2, fraction2
    """
    parsed_data: Dict = {}

    for line in f:
        line = line.strip()
        if not line:
            continue

        # Skip header line
        if line.startswith("#Quality"):
            continue

        parts = line.split("\t")
        if len(parts) < 5:
            continue

        try:
            quality = int(parts[0].lstrip("#"))
            count1 = int(parts[1])
            fraction1 = float(parts[2])
            count2 = int(parts[3])
            fraction2 = float(parts[4])
            parsed_data[quality] = [count1, fraction1, count2, fraction2]
        except (ValueError, IndexError):
            continue

    return parsed_data


def _prepare_plot_data(data_by_sample: Dict) -> list:
    """Prepare data for line graph plotting with count and proportion datasets."""
    # Calculate xmax based on 99.9% of data
    sumy = sum(
        int(data_by_sample[sample][x][0]) for sample in data_by_sample for x in data_by_sample[sample]
    )

    cutoff = sumy * 0.999
    all_x = set()
    for item in sorted(chain(*[data_by_sample[sample].items() for sample in data_by_sample])):
        all_x.add(item[0])
        cutoff -= item[1][0]
        if cutoff < 0:
            break

    # Columns to plot: Counts (0=Read1, 2=Read2), Proportions (1=Read1, 3=Read2)
    columns_to_plot = {
        "Counts": {0: "Read1", 2: "Read2"},
        "Proportions": {1: "Read1", 3: "Read2"},
    }

    plot_data = []
    for column_type in columns_to_plot:
        dataset = {}
        for column, column_name in columns_to_plot[column_type].items():
            for sample in data_by_sample:
                sample_key = f"{sample}.{column_name}"
                dataset[sample_key] = {}
                for x in all_x:
                    if x in data_by_sample[sample]:
                        dataset[sample_key][x] = data_by_sample[sample][x][column]
                    else:
                        dataset[sample_key][x] = 0
        plot_data.append(dataset)

    return plot_data
