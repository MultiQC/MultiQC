from itertools import chain

from multiqc.plots import linegraph

def plot_bhist(samples, file_type, **plot_args):
    """ Create line graph plot of histogram data for BBMap 'bhist' output.

    The 'samples' parameter could be from the bbmap mod_data dictionary:
    samples = bbmap.MultiqcModule.mod_data[file_type]
    """

    all_x = set()
    for item in sorted(chain(*[samples[sample]['data'].items()
                                for sample in samples])):
        all_x.add(item[0])

    columns_to_plot = {
        'GC': {
            1: 'C',
            2: 'G',
        },
        'AT': {
            0: 'A',
            3: 'T',
        },
        'N': {
            4: 'N'
        },
    }
    nucleotide_data = []
    for column_type in columns_to_plot:
        nucleotide_data.append(
            {
                sample+'.'+column_name: {
                    x: samples[sample]['data'][x][column]*100 if x in samples[sample]['data'] else 0
                    for x in all_x
            }
            for sample in samples
            for column, column_name in columns_to_plot[column_type].items()
        }
    )

    plot_params = {
            'id': 'bbmap-' + file_type + '_plot',
            'title': 'BBTools: ' + plot_args['plot_title'],
            'xlab': 'Read position',
            'ymin': 0,
            'ymax': 100,
            'data_labels': [
                {'name': 'Percentage of G+C bases'},
                {'name': 'Percentage of A+T bases'},
                {'name': 'Percentage of N bases'},
            ]
    }
    plot_params.update(plot_args['plot_params'])
    plot = linegraph.plot(
        nucleotide_data,
        plot_params
    )

    return plot
