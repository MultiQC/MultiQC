"""MultiQC submodule to parse BBTools indelhist output (indel length histogram)"""

import logging
from itertools import chain
from typing import Dict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


def parse_bbtools_indelhist(module: BaseMultiqcModule) -> int:
    """
    Parse BBTools indelhist output files showing indel length histogram.

    Shows the number of observed insertions and deletions for each length.
    """

    data_by_sample: Dict[str, Dict] = {}

    for f in module.find_log_files("bbtools/indelhist", filehandles=True):
        parsed = _parse_indelhist_file(f["f"])
        if parsed:
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="indelhist")
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
        "id": "bbtools-indelhist-plot",
        "title": "BBTools: Indel Lengths",
        "xlab": "Indel size",
        "xsuffix": " bp",
        "ylab": "Insertion count",
        "data_labels": [
            {"name": "Deletions", "ylab": "Deletion count"},
            {"name": "Insertions", "ylab": "Insertion count"},
        ],
    }

    module.add_section(
        name="Indel Lengths",
        anchor="bbtools-indelhist",
        description="Indel length histogram (`indelhist`). "
        "The plots show the number of observed insertions and deletions, "
        "for each insertion and deletion length.",
        plot=linegraph.plot(plot_data, pconfig),
    )

    # Write parsed data to file
    module.write_data_file(data_by_sample, "multiqc_bbtools_indelhist")

    return len(data_by_sample)


def _parse_indelhist_file(f) -> Dict:
    """
    Parse a BBTools indelhist file.

    Expected columns: #Length, Deletions, Insertions
    """
    parsed_data: Dict = {}

    for line in f:
        line = line.strip()
        if not line:
            continue

        # Skip header line
        if line.startswith("#Length"):
            continue

        parts = line.split("\t")
        if len(parts) < 3:
            continue

        try:
            length = int(parts[0].lstrip("#"))
            deletions = int(parts[1])
            insertions = int(parts[2])
            parsed_data[length] = [deletions, insertions]
        except (ValueError, IndexError):
            continue

    return parsed_data


def _prepare_plot_data(data_by_sample: Dict) -> list:
    """Prepare data for line graph plotting with Deletions and Insertions datasets."""
    # Get all x values
    all_x = set()
    for sample_data in data_by_sample.values():
        all_x.update(sample_data.keys())

    # Columns: Deletions = 0, Insertions = 1
    columns_to_plot = {
        "Deletions": {0: "Count"},
        "Insertions": {1: "Count"},
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
