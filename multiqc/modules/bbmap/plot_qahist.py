from itertools import chain

from multiqc.plots import linegraph


def plot_qahist(samples, file_type, **plot_args):
    """Create line graph plot of histogram data for BBMap 'qahist' output.

    The 'samples' parameter could be from the bbmap mod_data dictionary:
    samples = bbmap.MultiqcModule.mod_data[file_type]
    """

    all_x = set()
    for item in sorted(chain(*[samples[sample]["data"].items() for sample in samples])):
        # Skip plotting values for this x if missing columns.
        if len(item[1]) < 6:
            continue
        all_x.add(item[0])

    columns_to_plot = {
        "Match": {
            0: "Match",
        },
        "Sub": {
            1: "Sub",
        },
        "Ins": {
            2: "Ins",
        },
        "Del": {
            3: "Del",
        },
        "TrueQuality": {
            4: "TrueQuality",
        },
        "TrueQualitySub": {
            5: "TrueQualitySub",
        },
    }

    plot_data = []
    for column_type in columns_to_plot:
        column_data = {}
        for column, column_name in columns_to_plot[column_type].items():
            for sample in samples:
                sample_dict = {}
                for x in all_x:
                    if x in samples[sample]["data"]:
                        try:
                            sample_dict[x] = samples[sample]["data"][x][column]
                        except IndexError:
                            # This can happen if the data is missing columns.
                            # Which happens if e.g. a certain quality value is not present in the data.
                            sample_dict[x] = None
                    else:
                        sample_dict[x] = 0
                column_data[sample + "." + column_name] = sample_dict

        plot_data.append(column_data)
    plot_params = {
        "id": "bbmap-" + file_type + "_plot",
        "title": "BBTools: " + plot_args["plot_title"],
        "xlab": "Quality score",
        "ylab": "Match count",
        "data_labels": [
            {"name": "Match", "ylab": "Match count"},
            {"name": "Substitution", "ylab": "Substitution count"},
            {"name": "Insertion", "ylab": "Insertion count"},
            {"name": "Deletion", "ylab": "Deletion count"},
            {"name": "TrueQuality", "ylab": "Count"},
            {"name": "TrueQualitySubtitution", "ylab": "Count"},
        ],
    }

    plot_params.update(plot_args["plot_params"])
    plot = linegraph.plot(plot_data, plot_params)

    return plot
