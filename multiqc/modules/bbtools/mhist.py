"""MultiQC submodule to parse BBTools mhist output (match/sub/del/ins rates)"""

import logging
from itertools import chain
from typing import Dict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


def parse_bbtools_mhist(module: BaseMultiqcModule) -> int:
    """
    Parse BBTools mhist output files showing match, substitution, deletion, and insertion rates.

    Shows rates by read location.
    """

    data_by_sample: Dict[str, Dict] = {}

    for f in module.find_log_files("bbtools/mhist", filehandles=True):
        parsed = _parse_mhist_file(f["f"])
        if parsed:
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="mhist")
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
        "id": "bbtools-mhist-plot",
        "title": "BBTools: Match/Substitution/Deletion/Insertion Rates",
        "xlab": "Location in read",
        "xsuffix": " bp",
        "ylab": "Proportion",
        "data_labels": [
            {"name": "Read 1", "ylab": "Proportion"},
            {"name": "Read 2", "ylab": "Proportion"},
        ],
    }

    module.add_section(
        name="BBMap Match/Substitution/Deletion/Insertion Rates",
        anchor="bbtools-mhist",
        description="Histogram of match, substitution, deletion, and insertion rates by read location (`mhist`).",
        plot=linegraph.plot(plot_data, pconfig),
    )

    # Write parsed data to file
    module.write_data_file(data_by_sample, "multiqc_bbtools_mhist")

    return len(data_by_sample)


def _parse_mhist_file(f) -> Dict:
    """
    Parse a BBTools mhist file.

    Expected columns: #BaseNum, Match1, Sub1, Del1, Ins1, N1, Other1,
                      Match2, Sub2, Del2, Ins2, N2, Other2
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
        if len(parts) < 13:
            continue

        try:
            base_num = int(parts[0].lstrip("#"))
            match1 = float(parts[1])
            sub1 = float(parts[2])
            del1 = float(parts[3])
            ins1 = float(parts[4])
            n1 = float(parts[5])
            other1 = float(parts[6])
            match2 = float(parts[7])
            sub2 = float(parts[8])
            del2 = float(parts[9])
            ins2 = float(parts[10])
            n2 = float(parts[11])
            other2 = float(parts[12])
            parsed_data[base_num] = [
                match1,
                sub1,
                del1,
                ins1,
                n1,
                other1,
                match2,
                sub2,
                del2,
                ins2,
                n2,
                other2,
            ]
        except (ValueError, IndexError):
            continue

    return parsed_data


def _prepare_plot_data(data_by_sample: Dict) -> list:
    """Prepare data for line graph plotting with Read 1 and Read 2 datasets."""
    # Get all x values
    all_x = set()
    for sample_data in data_by_sample.values():
        all_x.update(sample_data.keys())

    # Read 1 columns: Match=0, Sub=1, Del=2, Ins=3, N1=4, Other1=5
    # Read 2 columns: Match=6, Sub=7, Del=8, Ins=9, N2=10, Other2=11
    columns_to_plot = {
        "Read 1": {
            0: "Match",
            1: "Sub",
            2: "Del",
            3: "Ins",
            4: "N1",
            5: "Other1",
        },
        "Read 2": {
            6: "Match",
            7: "Sub",
            8: "Del",
            9: "Ins",
            10: "N2",
            11: "Other2",
        },
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
