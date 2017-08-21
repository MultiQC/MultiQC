from itertools import chain

from multiqc.plots import linegraph

def plot_aqhist(samples, file_type, **plot_args):
    """ Create line graph plot of histogram data for BBMap 'aqhist' output.

    The 'samples' parameter could be from the bbmap mod_data dictionary:
    samples = bbmap.MultiqcModule.mod_data[file_type]
    """

    print(samples)

    sumy = sum([int(samples[sample]['data'][x][0])
                for sample in samples
                for x in samples[sample]['data']])

    cutoff = sumy * 0.999
    all_x = set()
    for item in sorted(chain(*[samples[sample]['data'].items()
                                for sample in samples])):
        print(item)
        all_x.add(item[0])
        cutoff -= item[1][0]
        if cutoff < 0:
            xmax = item[0]
            break

    columns_to_plot = [0,2]
    data = {
        sample+str(column): {
            x: samples[sample]['data'][x][column] if x in samples[sample]['data'] else 0
            for x in all_x
        }
        for sample in samples
        for column in columns_to_plot
    }

    plot_params = {
            'id': 'bbmap-' + file_type,
            'title': plot_args['plot_title'],
            'xmax': xmax
    }
    plot_params.update(plot_args['plot_params'])
    plot = linegraph.plot(
        data,
        plot_params
    )

    return plot