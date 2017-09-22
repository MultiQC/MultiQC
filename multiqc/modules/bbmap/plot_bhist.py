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
        'Nucleotides': {
            0: 'A',
            1: 'C',
            2: 'G',
            3: 'T',
            4: 'N'
        },
    }
    nucleotide_data = {
        sample+'.'+column_name: {
            x: samples[sample]['data'][x][column] if x in samples[sample]['data'] else 0
            for x in all_x
        }
        for sample in samples
        for column, column_name in columns_to_plot['Nucleotides'].items()
    }

    plot_params = {
            'id': 'bbmap-' + file_type,
            'title': plot_args['plot_title'],
            'xlab': 'Read position',
            'ylab': 'Proportion',
    }
    plot_params.update(plot_args['plot_params'])
    plot = linegraph.plot(
        nucleotide_data,
        plot_params
    )

    return plot