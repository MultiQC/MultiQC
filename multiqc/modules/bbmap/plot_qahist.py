from itertools import chain

from multiqc.plots import linegraph

def plot_qahist(samples, file_type, **plot_args):
    """ Create line graph plot of histogram data for BBMap 'qahist' output.

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
        # Skip plotting values for this x if missing columns.
        if len(item[1]) < 6:
            continue
        all_x.add(item[0])


    columns_to_plot = {
        'Match': {
            0: 'Match',
        },
        'Sub': {
            1: 'Sub',
        },
        'Ins': {
            2: 'Ins',
        },
        'Del': {
            3: 'Del',
        },
        'TrueQuality': {
            4: 'TrueQuality',
        },
        'TrueQualitySub': {
            5: 'TrueQualitySub',
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
            'xlab': 'Quality score',
            'data_labels': [
                {'name': 'Match', 'ylab': 'Match count'},
                {'name': 'Substitution', 'ylab': 'Substitution count'},
                {'name': 'Insertion', 'ylab': 'Insertion count'},
                {'name': 'Deletion', 'ylab': 'Deletion count'},
                {'name': 'TrueQuality', 'ylab': 'Count'},
                {'name': 'TrueQualitySubtitution', 'ylab': 'Count'},
            ]
    }

    plot_params.update(plot_args['plot_params'])
    plot = linegraph.plot(
        plot_data,
        plot_params
    )

    return plot
