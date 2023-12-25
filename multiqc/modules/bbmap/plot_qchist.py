from itertools import chain

from multiqc.plots import linegraph


def plot_qchist(samples, file_type, **plot_args):
    """Create line graph plot of histogram data for BBMap 'qchist' output.

    The 'samples' parameter could be from the bbmap mod_data dictionary:
    samples = bbmap.MultiqcModule.mod_data[file_type]
    """

    sumy = sum([int(samples[sample]["data"][x][0]) for sample in samples for x in samples[sample]["data"]])

    cutoff = sumy * 0.999
    all_x = set()
    for item in sorted(chain(*[samples[sample]["data"].items() for sample in samples])):
        all_x.add(item[0])
        cutoff -= item[1][0]
        if cutoff < 0:
            xmax = item[0]
            break
    else:
        xmax = max(all_x)

    data = {
        sample: {x: samples[sample]["data"][x][0] if x in samples[sample]["data"] else 0 for x in all_x}
        for sample in samples
    }
    # Add a count of 0.1 to zero counts, to avoid broken series in log axis
    data = {s: {k: d + 0.1 if d == 0 else d for k, d in v.items()} for s, v in data.items()}

    plot_params = {
        "id": "bbmap-" + file_type + "_plot",
        "title": "BBTools: " + plot_args["plot_title"],
        "xmax": xmax,
    }
    plot_params.update(plot_args["plot_params"])
    plot = linegraph.plot(data, plot_params)

    return plot
