from itertools import chain

from multiqc.plots import bargraph


def plot_bargraph(samples, file_type, **plot_args):
    """Create bar graph plot for BBDuk 'stats' output.

    The 'samples' parameter could be from the bbtools mod_data dictionary:
    samples = bbtools.MultiqcModule.mod_data[file_type]
    """

    all_x = set()
    for item in sorted(chain(*[samples[sample]["data"].items() for sample in samples])):
        all_x.add(item[0])

    columns_to_plot = {
        "Reads": {
            0: "Count",
        },
        "ReadsPct": {
            1: "Pct",
        },
    }

    plot_data = []
    for column_type in columns_to_plot:
        plot_data.append(
            {
                sample: {
                    x: (
                        samples[sample]["data"][x][column]
                        if x in samples[sample]["data"]
                        else 0
                    )
                    for x in all_x
                }
                for sample in samples
                for column, column_name in columns_to_plot[column_type].items()
            }
        )

    plot_params = {
        "id": f"bbtools-{file_type}_plot",
        "title": "BBTools: " + plot_args["plot_title"],
        "data_labels": [
            {"name": "Read Counts", "ylab": "Read Count", "tt_decimals": 0},
            {"name": "Read %", "ylab": "Read Percentage", "tt_decimals": 5},
        ],
    }
    plot_params.update(plot_args["plot_params"])
    plot = bargraph.plot(plot_data, pconfig=plot_params)

    return plot
