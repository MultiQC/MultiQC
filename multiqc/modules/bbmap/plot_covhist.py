from itertools import chain

from multiqc.plots import linegraph

def plot_covhist(samples, file_type, **plot_args):
    """ Create line graph plot for basic histogram data for 'covhist'.

    The 'samples' parameter could be from the bbmap mod_data dictionary:
    samples = bbmap.MultiqcModule.mod_data[file_type]
    """

    sumy = sum([int(samples[sample]['data'][x][0])
                for sample in samples
                for x in samples[sample]['data']])

    cutoff = sumy * 0.999
    all_x = set()
    for item in sorted(chain(*[samples[sample]['data'].items()
                                for sample in samples])):
        all_x.add(item[0])
        cutoff -= item[1][0]
        if cutoff < 0:
            xmax = item[0]
            break

    data = {
        sample: {
            x: samples[sample]['data'][x][0] if x in samples[sample]['data'] else 0
            for x in all_x
        }
        for sample in samples
    }

    plot_params = {
            'id': 'bbmap-' + file_type + '_plot',
            'title': 'BBTools: ' + plot_args['plot_title'],
            'smooth_points': 400,
            'xmax': xmax,
            'xlab': 'Coverage (depth)',
            'ylab': 'Number of occurences'
    }
    plot_params.update(plot_args['plot_params'])
    plot = linegraph.plot(
        data,
        plot_params
    )

    return plot
