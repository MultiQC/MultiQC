"""MultiQC submodule to parse BBTools bqhist output (base quality boxplot data)"""

import logging
from itertools import chain
from typing import Dict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


def parse_bbtools_bqhist(module: BaseMultiqcModule) -> int:
    """
    Parse BBTools bqhist output files - quality histogram designed for box plots.

    Shows mean base quality for each read position.
    """

    data_by_sample: Dict[str, Dict] = {}

    for f in module.find_log_files("bbtools/bqhist", filehandles=True):
        parsed = _parse_bqhist_file(f["f"])
        if parsed:
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="bqhist")
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
        "id": "bbtools-bqhist-plot",
        "title": "BBTools: Base Quality",
        "xlab": "Read position",
        "ylab": "Average quality score",
        "data_labels": [
            {"name": "Read 1"},
            {"name": "Read 2"},
        ],
    }

    module.add_section(
        name="BBMap Base Quality",
        anchor="bbtools-bqhist",
        description="Quality histogram designed for box plots (`bqhist`). "
        "Refer to original source files for complete boxplot data. "
        "Plot shows mean base quality for each read position.",
        plot=linegraph.plot(plot_data, pconfig),
    )

    # Write parsed data to file
    module.write_data_file(data_by_sample, "multiqc_bbtools_bqhist")

    return len(data_by_sample)


def _parse_bqhist_file(f) -> Dict:
    """
    Parse a BBTools bqhist file.

    Expected columns: #BaseNum, count_1, min_1, max_1, mean_1, Q1_1, med_1, Q3_1, LW_1, RW_1,
                      count_2, min_2, max_2, mean_2, Q1_2, med_2, Q3_2, LW_2, RW_2
    """
    parsed_data: Dict = {}

    for line in f:
        line = line.strip()
        if not line:
            continue

        # Skip header line
        if line.startswith("#BaseNum"):
            continue

        parts = line.split("\t")
        if len(parts) < 18:
            continue

        try:
            base_num = int(parts[0].lstrip("#"))
            count_1 = int(parts[1])
            min_1 = int(parts[2])
            max_1 = int(parts[3])
            mean_1 = float(parts[4])
            q1_1 = int(parts[5])
            med_1 = int(parts[6])
            q3_1 = int(parts[7])
            lw_1 = int(parts[8])
            rw_1 = int(parts[9])
            count_2 = int(parts[10])
            min_2 = int(parts[11])
            max_2 = int(parts[12])
            mean_2 = float(parts[13])
            q1_2 = int(parts[14])
            med_2 = int(parts[15])
            q3_2 = int(parts[16])
            lw_2 = int(parts[17])
            rw_2 = int(parts[18]) if len(parts) > 18 else 0
            parsed_data[base_num] = [
                count_1, min_1, max_1, mean_1, q1_1, med_1, q3_1, lw_1, rw_1,
                count_2, min_2, max_2, mean_2, q1_2, med_2, q3_2, lw_2, rw_2,
            ]
        except (ValueError, IndexError):
            continue

    return parsed_data


def _prepare_plot_data(data_by_sample: Dict) -> list:
    """Prepare data for line graph plotting with Read 1 and Read 2 mean quality."""
    # Get all x values
    all_x = set()
    for sample_data in data_by_sample.values():
        all_x.update(sample_data.keys())

    # Mean quality columns: Read 1 = index 3, Read 2 = index 12
    columns_to_plot = {
        "Read 1 averages": {3: "Mean"},
        "Read 2 averages": {12: "Mean"},
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
