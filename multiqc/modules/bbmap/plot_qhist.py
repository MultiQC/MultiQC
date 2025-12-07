from itertools import chain
from typing import Any, Dict, List, Set

from multiqc.plots import linegraph


def plot_qhist(samples: Dict[str, Any], file_type: str, **plot_args: Any):
    """Create line graph plot of histogram data for BBMap 'qhist' output.

    The 'samples' parameter could be from the bbmap mod_data dictionary:
    samples = bbmap.MultiqcModule.mod_data[file_type]
    """

    sumy = sum([int(samples[sample]["data"][x][0]) for sample in samples for x in samples[sample]["data"]])

    cutoff = sumy * 0.999
    all_x: Set[int] = set()
    for item in sorted(chain(*[samples[sample]["data"].items() for sample in samples])):
        all_x.add(item[0])
        cutoff -= item[1][0]
        if cutoff < 0:
            xmax = item[0]
            break
    else:
        xmax = max(all_x)

    columns_to_plot: Dict[str, Dict[int, str]] = {
        "Linear": {0: "Read1", 3: "Read2"},
        "Logarithmic": {1: "Read1", 4: "Read2"},
        "Measured": {2: "Read1", 5: "Read2"},
    }

    plot_data: List[Dict[str, Any]] = []
    for column_type in columns_to_plot:
        plot_data.append(
            {
                sample + "." + column_name: {
                    x: samples[sample]["data"][x][column] if x in samples[sample]["data"] else 0
                    for x in all_x
                    if len(samples[sample]["data"]) > x and len(samples[sample]["data"][x]) > column
                }
                for sample in samples
                for column, column_name in columns_to_plot[column_type].items()
            }
        )

    plot_params = {
        "id": "bbmap-" + file_type + "_plot",
        "title": "BBTools: " + plot_args["plot_title"],
        "xmax": xmax,
        "xlab": "Position in read",
        "xsuffix": " bp",
        "ylab": "Quality score",
        "data_labels": [
            {"name": "Linear", "ylab": "Quality score"},
            {"name": "Logarithmic", "ylab": "Log score"},
            {"name": "Measured", "ylab": "Measured quality value"},
        ],
    }

    plot_params.update(plot_args["plot_params"])
    plot = linegraph.plot(plot_data, plot_params)

    return plot
