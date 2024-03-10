from itertools import chain
from typing import Dict

from multiqc.plots import linegraph


def plot_qahist(data_by_sample: Dict[str, Dict[str, Dict]], file_type, **plot_args):
    """
    Create line graph plot of histogram data for BBMap 'qahist' output.

    The 'samples' parameter could be from the bbmap mod_data dictionary:
    samples = bbmap.MultiqcModule.mod_data[file_type]

    samples:
        sample1:
            data:
                0: [128, 0, 0, 60]
                1: ...
            kv:
                Deviation: '13.882'
                DeviationSub: '13.056'
                ...
        sample2:
            ...
    """

    all_x = set()
    for item in sorted(chain(*[sample_data["data"].items() for sample, sample_data in data_by_sample.items()])):
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
        y_by_x_by_sample = {}
        for column, column_name in columns_to_plot[column_type].items():
            for sample, sample_data in data_by_sample.items():
                y_by_x = {}
                for x in all_x:
                    if x not in sample_data["data"] or column >= len(sample_data["data"][x]):
                        # This can happen if the data is missing columns.
                        # Which happens if e.g. a certain quality value is not present in the data.
                        pass
                    else:
                        y_by_x[x] = sample_data["data"][x][column]
                y_by_x_by_sample[sample] = y_by_x

        plot_data.append(y_by_x_by_sample)
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
            {"name": "TrueQualitySubstitution", "ylab": "Count"},
        ],
        "xmin": 0,
        "ymin": 0,
    }

    plot_params.update(plot_args["plot_params"])
    plot = linegraph.plot(plot_data, plot_params)

    return plot
