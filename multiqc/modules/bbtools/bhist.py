"""MultiQC submodule to parse BBTools bhist output (base composition histogram)"""

import logging
from itertools import chain
from typing import Dict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


def parse_bbtools_bhist(module: BaseMultiqcModule) -> int:
    """
    Parse BBTools bhist output files showing base composition histogram by position.

    Shows the percentage of G+C, A+T, and N bases for each position in the reads.
    """

    data_by_sample: Dict[str, Dict] = {}

    for f in module.find_log_files("bbtools/bhist", filehandles=True):
        parsed = _parse_bhist_file(f["f"])
        if parsed:
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="bhist")
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
        "id": "bbtools-bhist-plot",
        "title": "BBTools: Base Composition",
        "xlab": "Read position",
        "xsuffix": " bp",
        "ylab": "Percentage of G+C bases",
        "ymin": 0,
        "ymax": 100,
        "data_labels": [
            {"name": "Percentage of G+C bases"},
            {"name": "Percentage of A+T bases"},
            {"name": "Percentage of N bases"},
        ],
    }

    module.add_section(
        name="Base Composition",
        anchor="bbtools-bhist",
        description="Base composition histogram by position (`bhist`). "
        "The plot shows the percentage of `G+C`, `A+T`, and `N` bases "
        "for each position in the reads.",
        helptext="Relative composition",
        plot=linegraph.plot(plot_data, pconfig),
    )

    # Write parsed data to file
    module.write_data_file(data_by_sample, "multiqc_bbtools_bhist")

    return len(data_by_sample)


def _parse_bhist_file(f) -> Dict:
    """
    Parse a BBTools bhist file.

    Expected columns: #Pos, A, C, G, T, N
    """
    parsed_data: Dict = {}

    for line in f:
        line = line.strip()
        if not line:
            continue

        # Skip header line
        if line.startswith("#Pos"):
            continue

        parts = line.split("\t")
        if len(parts) < 6:
            continue

        try:
            pos = int(parts[0].lstrip("#"))
            a = float(parts[1])
            c = float(parts[2])
            g = float(parts[3])
            t = float(parts[4])
            n = float(parts[5])
            parsed_data[pos] = [a, c, g, t, n]
        except (ValueError, IndexError):
            continue

    return parsed_data


def _prepare_plot_data(data_by_sample: Dict) -> list:
    """Prepare data for line graph plotting with GC, AT, and N datasets."""
    # Get all x values
    all_x = set()
    for sample_data in data_by_sample.values():
        all_x.update(sample_data.keys())

    # Columns: A=0, C=1, G=2, T=3, N=4
    columns_to_plot = {
        "GC": {1: "C", 2: "G"},
        "AT": {0: "A", 3: "T"},
        "N": {4: "N"},
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
                        # Multiply by 100 to get percentage
                        dataset[sample_key][x] = data_by_sample[sample][x][column] * 100
                    else:
                        dataset[sample_key][x] = 0
        plot_data.append(dataset)

    return plot_data
