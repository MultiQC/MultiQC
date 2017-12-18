from itertools import chain

from multiqc.plots import linegraph

def plot_indelhist(samples, file_type, **plot_args):
    """ Create line graph plot of histogram data for BBMap 'indelhist' output.

    The 'samples' parameter could be from the bbmap mod_data dictionary:
    samples = bbmap.MultiqcModule.mod_data[file_type]
    """

    all_x = set()
    for item in sorted(chain(*[samples[sample]['data'].items()
                                for sample in samples])):
        all_x.add(item[0])


    columns_to_plot = {
        'Deletions': {
            0: 'Count',
        },
        'Insertions': {
            1: 'Count',
        }
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
            'xlab': 'Indel size',
            'data_labels': [
                {'name': 'Insertions', 'ylab': 'Insertion count'},
                {'name': 'Deletions', 'ylab': 'Deletion count'},
            ]

    }
    plot_params.update(plot_args['plot_params'])
    plot = linegraph.plot(
        plot_data,
        plot_params
    )

    return plot
