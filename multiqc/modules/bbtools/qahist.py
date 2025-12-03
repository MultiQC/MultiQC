"""MultiQC submodule to parse BBTools qahist output (quality accuracy histogram)"""

import logging
from itertools import chain
from typing import Dict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, table

log = logging.getLogger(__name__)


def parse_bbtools_qahist(module: BaseMultiqcModule) -> int:
    """
    Parse BBTools qahist output files showing quality accuracy histogram.

    Shows observed count of each type of alignment by base quality score.
    """

    data_by_sample: Dict[str, Dict] = {}

    for f in module.find_log_files("bbtools/qahist", filehandles=True):
        parsed = _parse_qahist_file(f["f"])
        if parsed:
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="qahist")
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
        "id": "bbtools-qahist-plot",
        "title": "BBTools: Quality Accuracy",
        "xlab": "Quality score",
        "ylab": "Match count",
        "xmin": 0,
        "ymin": 0,
        "data_labels": [
            {"name": "Match", "ylab": "Match count"},
            {"name": "Substitution", "ylab": "Substitution count"},
            {"name": "Insertion", "ylab": "Insertion count"},
            {"name": "Deletion", "ylab": "Deletion count"},
            {"name": "TrueQuality", "ylab": "Count"},
            {"name": "TrueQualitySubstitution", "ylab": "Count"},
        ],
    }

    module.add_section(
        name="Quality Accuracy",
        anchor="bbtools-qahist",
        description="Base quality accuracy histogram of error rates versus quality score (`qahist`). "
        "The plots show the observed count of each type of alignment by base quality score.",
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
            "DeviationSub": {
                "title": "Deviation (Sub)",
                "description": "Deviation from expected quality (substitutions)",
            },
        }

        module.add_section(
            name="Quality Accuracy Summary",
            anchor="bbtools-qahist-table",
            description="Quality accuracy statistics summary.",
            plot=table.plot(
                table_data,
                table_headers,
                pconfig={
                    "id": "bbtools-qahist-summary-table",
                    "namespace": "BBTools",
                    "title": "BBTools: Quality Accuracy Summary",
                },
            ),
        )

    # Write parsed data to file
    module.write_data_file(data_by_sample, "multiqc_bbtools_qahist")

    return len(data_by_sample)


def _parse_qahist_file(f) -> Dict:
    """
    Parse a BBTools qahist file.

    Expected columns: #Quality, Match, Sub, Ins, Del, TrueQuality, TrueQualitySub
    Also contains key-value pairs like Deviation, DeviationSub
    """
    parsed_data: Dict = {"data": {}, "kv": {}}
    kv_keys = ["Deviation", "DeviationSub"]

    for line in f:
        line = line.strip()
        if not line:
            continue

        parts = line.split("\t")

        # Check for header/kv lines
        if parts[0].startswith("#"):
            parts[0] = parts[0][1:]  # Remove leading '#'

            # Skip header
            if parts[0] == "Quality":
                continue

            # Key-value pair
            if len(parts) == 2 and parts[0] in kv_keys:
                try:
                    parsed_data["kv"][parts[0]] = float(parts[1])
                except ValueError:
                    parsed_data["kv"][parts[0]] = parts[1]
                continue

        # Data row
        if len(parts) >= 6:
            try:
                quality = int(parts[0].lstrip("#"))
                match_val = int(parts[1])
                sub = int(parts[2])
                ins = int(parts[3])
                delete = int(parts[4])
                true_quality = float(parts[5]) if len(parts) > 5 else 0.0
                true_quality_sub = float(parts[6]) if len(parts) > 6 else 0.0
                parsed_data["data"][quality] = [match_val, sub, ins, delete, true_quality, true_quality_sub]
            except (ValueError, IndexError):
                continue

    return parsed_data if parsed_data["data"] else {}


def _prepare_plot_data(data_by_sample: Dict) -> list:
    """Prepare data for line graph plotting with multiple data series."""
    # Get all x values, but only include those with complete data
    all_x = set()
    for sample_data in data_by_sample.values():
        for x, vals in sample_data.get("data", {}).items():
            if len(vals) >= 6:
                all_x.add(x)

    # Columns: Match=0, Sub=1, Ins=2, Del=3, TrueQuality=4, TrueQualitySub=5
    columns_to_plot = {
        "Match": {0: "Match"},
        "Sub": {1: "Sub"},
        "Ins": {2: "Ins"},
        "Del": {3: "Del"},
        "TrueQuality": {4: "TrueQuality"},
        "TrueQualitySub": {5: "TrueQualitySub"},
    }

    plot_data = []
    for column_type in columns_to_plot:
        dataset = {}
        for column, column_name in columns_to_plot[column_type].items():
            for sample in data_by_sample:
                dataset[sample] = {}
                sample_data = data_by_sample[sample].get("data", {})
                for x in all_x:
                    if x in sample_data and column < len(sample_data[x]):
                        dataset[sample][x] = sample_data[x][column]
                    # Don't add missing values - leave gaps in plot
        plot_data.append(dataset)

    return plot_data
