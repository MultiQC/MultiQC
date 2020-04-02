from itertools import chain

from multiqc.plots import linegraph

def plot_mhist(samples, file_type, **plot_args):
    """ Create line graph plot of histogram data for BBMap 'mhist' output.

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


    columns_to_plot = {
        'Read 1': {
            0: 'Match',
            1: 'Sub',
            2: 'Del',
            3: 'Ins',
            4: 'N1',
            5: 'Other1',
        },
        'Read 2': {
            6: 'Match',
            7: 'Sub',
            8: 'Del',
            9: 'Ins',
            10: 'N2',
            11: 'Other2',
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
            'xlab': 'Location in read',
            'data_labels': [
                {'name': 'Read 1', 'ylab': 'Proportion'},
                {'name': 'Read 2', 'ylab': 'Proportion'},
            ]
    }

    plot_params.update(plot_args['plot_params'])
    plot = linegraph.plot(
        plot_data,
        plot_params
    )

    return plot
