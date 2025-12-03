"""MultiQC submodule to parse BBTools qchist output (quality value counts)"""

import logging
from itertools import chain
from typing import Dict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

log = logging.getLogger(__name__)

# Column definitions for qchist file
# Quality, count1, fraction1
QCHIST_COLS = {
    "Quality": int,
    "count1": int,
    "fraction1": float,
}


def parse_bbtools_qchist(module: BaseMultiqcModule) -> int:
    """
    Parse BBTools qchist output files showing count of bases with each quality value.

    Also calculates Q30 percentage and adds it to general stats.
    """

    data_by_sample: Dict[str, Dict] = {}

    for f in module.find_log_files("bbtools/qchist", filehandles=True):
        parsed = _parse_qchist_file(f["f"])
        if parsed:
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="qchist")
            data_by_sample[s_name] = parsed

    # Filter to strip out ignored sample names
    data_by_sample = module.ignore_samples(data_by_sample)

    if len(data_by_sample) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    module.add_software_version(None)

    # Calculate Q30 percentage for general stats
    general_stats_data: Dict[str, Dict] = {}
    for s_name, sample_data in data_by_sample.items():
        fraction_gt_q30 = []
        for qual, values in sample_data.items():
            if int(qual) >= 30:
                # fraction1 is at index 1 (after count1)
                fraction_gt_q30.append(values[1])
        general_stats_data[s_name] = {"pct_q30": sum(fraction_gt_q30) * 100.0}

    general_stats_headers: Dict = {
        "pct_q30": {
            "title": "% Q30 bases",
            "description": "BBTools qchist - Percentage of bases with phred quality score >= 30",
            "suffix": "%",
            "scale": "RdYlGn",
            "format": "{:,.2f}",
            "min": 0,
            "max": 100,
        }
    }

    general_stats_headers = module.get_general_stats_headers(all_headers=general_stats_headers)
    if general_stats_headers:
        module.general_stats_addcols(general_stats_data, general_stats_headers, namespace="qchist")

    # Create line graph plot
    plot_data = _prepare_plot_data(data_by_sample)

    pconfig = {
        "id": "bbtools-qchist-plot",
        "title": "BBTools: Count of Bases with Each Quality Value",
        "xlab": "Phred Score",
        "ylab": "Counts",
        "ylog": True,
        "x_bands": [
            {"from": 30, "to": 100, "color": "#009500", "opacity": 0.13},
            {"from": 20, "to": 30, "color": "#a07300", "opacity": 0.13},
            {"from": 0, "to": 20, "color": "#990101", "opacity": 0.13},
        ],
    }

    module.add_section(
        name="Count of Bases with Each Quality Value",
        anchor="bbtools-qchist",
        description="Histogram of base qualities (`qchist`). "
        "Plot shows the number of bases at each quality score. Zero counts are shown as `0.1` due to log axis.",
        plot=linegraph.plot(plot_data, pconfig),
    )

    # Write parsed data to file
    module.write_data_file(data_by_sample, "multiqc_bbtools_qchist")

    return len(data_by_sample)


def _parse_qchist_file(f) -> Dict:
    """
    Parse a BBTools qchist file.

    Expected columns: #Quality, count1, fraction1
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
        if len(parts) < 3:
            continue

        try:
            quality = int(parts[0].lstrip("#"))
            count1 = int(parts[1])
            fraction1 = float(parts[2])
            parsed_data[quality] = [count1, fraction1]
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

    # Build plot data (count1 values)
    plot_data = {}
    for sample in data_by_sample:
        plot_data[sample] = {}
        for x in all_x:
            if x in data_by_sample[sample]:
                count = data_by_sample[sample][x][0]
                # Add 0.1 to zero counts to avoid broken series in log axis
                plot_data[sample][x] = count + 0.1 if count == 0 else count
            else:
                plot_data[sample][x] = 0.1

    return plot_data
