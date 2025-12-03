"""MultiQC submodule to parse BBTools ehist output (errors-per-read histogram)"""

import logging
from itertools import chain
from typing import Dict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


def parse_bbtools_ehist(module: BaseMultiqcModule) -> int:
    """
    Parse BBTools ehist output files showing errors-per-read histogram.
    """

    data_by_sample: Dict[str, Dict] = {}

    for f in module.find_log_files("bbtools/ehist", filehandles=True):
        parsed = _parse_ehist_file(f["f"])
        if parsed:
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="ehist")
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
        "id": "bbtools-ehist-plot",
        "title": "BBTools: Errors-per-Read",
        "xlab": "Errors",
        "ylab": "# Reads",
    }

    module.add_section(
        name="BBMap Errors per Read",
        anchor="bbtools-ehist",
        description="Errors-per-read histogram (`ehist`).",
        plot=linegraph.plot(plot_data, pconfig),
    )

    # Write parsed data to file
    module.write_data_file(data_by_sample, "multiqc_bbtools_ehist")

    return len(data_by_sample)


def _parse_ehist_file(f) -> Dict:
    """
    Parse a BBTools ehist file.

    Expected columns: #Errors, Count
    """
    parsed_data: Dict = {}

    for line in f:
        line = line.strip()
        if not line:
            continue

        # Skip header line
        if line.startswith("#Errors"):
            continue

        parts = line.split("\t")
        if len(parts) < 2:
            continue

        try:
            errors = int(parts[0].lstrip("#"))
            count = int(parts[1])
            parsed_data[errors] = [count]
        except (ValueError, IndexError):
            continue

    return parsed_data


def _prepare_plot_data(data_by_sample: Dict) -> Dict:
    """Prepare data for line graph plotting."""
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

    # Build plot data
    plot_data = {}
    for sample in data_by_sample:
        plot_data[sample] = {}
        for x in all_x:
            if x in data_by_sample[sample]:
                plot_data[sample][x] = data_by_sample[sample][x][0]
            else:
                plot_data[sample][x] = 0

    return plot_data
