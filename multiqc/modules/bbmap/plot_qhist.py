from itertools import chain

from multiqc.plots import linegraph

def plot_qhist(samples, file_type, **plot_args):
    """ Create line graph plot of histogram data for BBMap 'qhist' output.

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

    # The columns_to_plot dictionary should be replaced with an OrderedDict,
    # not to rely on the ordered-by-default implementation details on Python 3.6
    columns_to_plot = {
        'Linear': {
            0: 'Read1',
            3: 'Read2'
        },
        'Logarithmic': {
            1: 'Read1',
            4: 'Read2'
        },
        'Measured': {
            2: 'Read1',
            5: 'Read2'
        },
    }

    plot_data = []
    for column_type in columns_to_plot:
        plot_data.append(
            {
                sample+'.'+column_name: {
                    x: samples[sample]['data'][x][column] if x in samples[sample]['data'] else 0
                    for x in all_x
                }
                for sample in samples
                for column, column_name in columns_to_plot[column_type].items()
            }
        )

    plot_params = {
            'id': 'bbmap-' + file_type + '_plot',
            'title': 'BBTools: ' + plot_args['plot_title'],
            'xmax': xmax,
            'xlab': 'Position in read',
            'data_labels': [
                {'name': 'Linear', 'ylab': 'Quality score'},
                {'name': 'Logarithmic', 'ylab': 'Log score'},
                {'name': 'Measured', 'ylab': 'Measured quality value'},
            ]
    }

    plot_params.update(plot_args['plot_params'])
    plot = linegraph.plot(
        plot_data,
        plot_params
    )

    return plot
