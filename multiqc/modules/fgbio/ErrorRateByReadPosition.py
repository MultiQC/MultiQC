#!/usr/bin/env python

""" MultiQC submodule to parse output from fgbio ErrorRateByReadPosition """


from multiqc.plots import linegraph


def parse_reports(self):
    """Parser metric files for fgbio ErrorRateByReadPosition.
    
    Stores the per-read-per-position metrics into a data file and adds a section
    with a per-sample plot.
    """

    # slurp in all the data
    all_data = dict()
    error_rates = dict()
    y_max = 0.01  # default to 1%
    for f in self.find_log_files('fgbio/errorratebyreadposition', filehandles=True):
        fh = f['f']
        header = fh.readline().rstrip('\r\n').split('\t')
        if not header or header[0] != 'read_number':
            continue
        
        # slurp in the data for this sample
        s_name = f['s_name']
        s_data = dict()
        bases_total = 0
        errors = 0
        for line in fh:
            fields = line.rstrip('\r\n').split('\t')
            fields[1:4] = [int(field) for field in fields[1:4]]
            fields[4:] = [float(field) for field in fields[4:]]
            row_data = dict(zip(header, fields))
            read_number = row_data['read_number']
            position = row_data['position']
            if read_number not in s_data:
                s_data[read_number] = dict()
            s_data[read_number][position] = row_data
            if row_data['error_rate'] > y_max:
                y_max = row_data['error_rate']
            bases_total += row_data['bases_total']
            errors += row_data['errors']

        if s_data:
            all_data[s_name] = s_data
            error_rate = 0.0 if bases_total == 0 else errors / float(bases_total)
            error_rates[s_name] = {'error_rate': error_rate}


    # ignore samples
    all_data = self.ignore_samples(all_data)

    # if no data, then do nothing
    if not all_data:
        return 0

    # Write parsed data to a file
    self.write_data_file(all_data, 'multiqc_fgbio_ErrorRateByReadPosition_per_position')
    self.write_data_file(error_rates, 'multiqc_fgbio_ErrorRateByReadPosition_total')

    # Plot the data and add section
    pconfig = {
        'id': 'fgbio_ErrorRateByReadPosition',
        'title': 'Fgbio: Error Rate by Read Position',
        'ylab': 'Error rate',
        'xlab': 'Read Position',
        'xDecimals': False,
        'tt_label': '<b>read position {point.x}</b>: {point.y:.2f} %',
        'ymax': y_max,
        'ymin': 0,
        'data_labels': [
            {'name': 'Error Rate', 'ylab': 'Error_Rate'},
            {'name': 'A to C error rate', 'ylab': 'A to C error rate'},
            {'name': 'A to G error rate', 'ylab': 'A to G error rate'},
            {'name': 'A to T error rate', 'ylab': 'A to T error rate'},
            {'name': 'C to A error rate', 'ylab': 'C to A error rate'},
            {'name': 'C to G error rate', 'ylab': 'C to G error rate'},
            {'name': 'C to T error rate', 'ylab': 'C to T error rate'},
        ]
    }

    # Build a list of linegraphs
    linegraph_data = [{}, {}, {}, {}, {}, {}, {}]
    linegraph_keys = ['error_rate', 'a_to_c_error_rate', 'a_to_g_error_rate', 'a_to_t_error_rate', 'c_to_a_error_rate', 'c_to_g_error_rate', 'c_to_t_error_rate']
    for s_name, s_data in all_data.items():
        for read_number, read_data in s_data.items():
            s_name = "%s_R%d" % (s_name, int(read_number) + 1)
            for lg, index in zip(linegraph_data, range(7)):
                lg[s_name] = dict((d['position'], d[linegraph_keys[index]]) for d in read_data.values())

    # add a section for the plot
    self.add_section (
        name = 'Error Rate by Read Position',
        anchor = 'fgbio-error-rate-by-read-position',
        description = 'Plot shows the error rate by read position.',
        plot = linegraph.plot(linegraph_data, pconfig)
    )

    # Add to general stats table
    y_max = 0.02
    for s_name in error_rates:
        if error_rates[s_name]['error_rate'] > y_max:
            y_max = error_rates[s_name]['error_rate']
    self.general_stats_headers['error_rate'] = {
        'title': 'Percent Error',
        'description': 'Percent error across all read positions',
        'min': 0,
        'max': y_max * 100.0,
        'scale': 'GnYlRd',
        'suffix': '%',
        'format': '{:,.2f}',
        'modify': lambda x: 100.0 * x
    }
    for s_name in error_rates:
        if s_name not in self.general_stats_data:
            self.general_stats_data[s_name] = dict()
        self.general_stats_data[s_name].update(error_rates[s_name])

    return len(all_data)
