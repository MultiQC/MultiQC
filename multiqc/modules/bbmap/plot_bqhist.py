from itertools import chain

from multiqc.plots import linegraph


def plot_bqhist(samples, file_type, **plot_args):
    """Create line graph plot of histogram data for BBMap 'bqhist' output.

    The 'samples' parameter could be from the bbmap mod_data dictionary:
    samples = bbmap.MultiqcModule.mod_data[file_type]
    """

    all_x = set()
    for item in sorted(chain(*[samples[sample]["data"].items() for sample in samples])):
        all_x.add(item[0])

    columns_to_plot = {
        #'Counts': {
        #    0: 'Read 1',
        #    9: 'Read 2',
        # },
        "Read 1 averages": {
            # 1: 'Min',
            # 2: 'Max',
            3: "Mean",
            # 4: '1st quartile',
            # 5: 'Median',
            # 6: '3rd quartile',
            # 7: 'Left (lower) whisker',
            # 8: 'Right (upper) whisker',
        },
        "Read 2 averages": {
            # 10: 'Min',
            # 11: 'Max',
            12: "Mean",
            # 13: '1st quartile',
            # 14: 'Median',
            # 15: '3rd quartile',
            # 16: 'Left (lower) whisker',
            # 17: 'Right (upper) whisker',
        },
    }

    plot_data = []
    for column_type in columns_to_plot:
        plot_data.append(
            {
                sample + "." + column_name: {
                    x: samples[sample]["data"][x][column] if x in samples[sample]["data"] else 0 for x in all_x
                }
                for sample in samples
                for column, column_name in columns_to_plot[column_type].items()
            }
        )

    plot_params = {
        "id": "bbmap-" + file_type + "_plot",
        "title": "BBTools: " + plot_args["plot_title"],
        "xlab": "Read position",
        "ylab": "Average quality score",
        "data_labels": [
            # {'name': 'Count histogram', 'ylab': 'Read count'},
            {
                "name": "Read 1",
            },
            {
                "name": "Read 2",
            },
        ],
    }

    plot_params.update(plot_args["plot_params"])
    plot = linegraph.plot(plot_data, plot_params)

    return plot
