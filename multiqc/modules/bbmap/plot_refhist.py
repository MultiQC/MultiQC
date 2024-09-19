from itertools import chain
from collections import OrderedDict

from multiqc.plots import bargraph


def plot_refhist(samples, file_type, **plot_args):
    """Create bargraph plot for basic histogram data for 'file_type'.

    The 'samples' parameter could be from the bbmap mod_data dictionary:
    samples = bbmap.MultiqcModule.mod_data[file_type]
    """

    all_x = set()
    for item in sorted(chain(*[samples[sample]["data"].items() for sample in samples])):
        all_x.add(item[0])

    columns_to_plot = {
        0: "Unambiguous Reads (%)",
        1: "Unambiguous Reads (MB)",
        4: "Unambiguous Reads",
        2: "Ambiguous Reads (%)",
        3: "Ambiguous Reads (MB)",
        5: "Ambiguous Reads",
        6: "Assigned Reads",
        7: "Assigned Bases",
    }

    plot_data = []
    for column, column_name in columns_to_plot.items():
        plot_data.append(
            {
                sample: {x: samples[sample]["data"][x][column] if x in samples[sample]["data"] else 0 for x in all_x}
                for sample in samples
            }
        )

    data_labels = [
        {"name": "Unambiguous Reads (%)", "ylab": "Percentage", "ymax": 100},
        {"name": "Unambiguous Reads (MB)", "ylab": "Megabases"},
        {"name": "Unambiguous Reads", "ylab": "Read count"},
        {"name": "Ambiguous Reads (%)", "ylab": "Percentage", "ymax": 100},
        {"name": "Ambiguous Reads (MB)", "ylab": "Megabases"},
        {"name": "Ambiguous Reads", "ylab": "Read count"},
        {"name": "Assigned Reads", "ylab": "Read count"},
        {"name": "Assigned Bases", "ylab": "Base count"},
    ]

    plot_params = {
        "id": "bbmap-" + file_type + "_plot",
        "title": "BBTools: " + plot_args["plot_title"],
        "xlab": "Sample",
        "ylab": "",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "data_labels": data_labels,
        "tt_decimals": 2,
        "tt_percentages": False,
        "stacking": None,
    }

    plot_params.update(plot_args["plot_params"])
    plot = bargraph.plot(plot_data, pconfig=plot_params)

    return plot
