#!/usr/bin/env python

""" MultiQC submodule to parse output from fgbio ErrorRateByReadPosition """


from multiqc.plots import linegraph


def parse_reports(self):
    """Parser metric files for fgbio ErrorRateByReadPosition.

    Stores the per-read-per-position metrics into a data file and adds a section
    with a per-sample plot.
    """
    linegraph_keys = ['error_rate', 'a_to_c_error_rate', 'a_to_g_error_rate', 'a_to_t_error_rate', 'c_to_a_error_rate', 'c_to_g_error_rate', 'c_to_t_error_rate']

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
            for key in linegraph_keys:
                y_max = max(y_max, row_data[key])
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
            {'name': 'Error Rate', 'ylab': 'Error Rate'},
            {'name': 'A > C', 'ylab': 'A to C error rate'},
            {'name': 'A > G', 'ylab': 'A to G error rate'},
            {'name': 'A > T', 'ylab': 'A to T error rate'},
            {'name': 'C > A', 'ylab': 'C to A error rate'},
            {'name': 'C > G', 'ylab': 'C to G error rate'},
            {'name': 'C > T', 'ylab': 'C to T error rate'},
        ]
    }

    # Build a list of linegraphs
    linegraph_data = [{} for _ in linegraph_keys]
    for s_name, s_data in all_data.items():
        for read_number, read_data in s_data.items():
            s_name_with_read = "%s_R%d" % (s_name, int(read_number))
            for lg, index in zip(linegraph_data, range(7)):
                lg[s_name_with_read] = dict((d['position'], d[linegraph_keys[index]]) for d in read_data.values())

    # add a section for the plot
    self.add_section (
        name = 'Error Rate by Read Position',
        anchor = 'fgbio-error-rate-by-read-position',
        description = 'Error rate by read position. Plot tabs show the error rates for specific substitution types.',
        helptext = '''
        The error rate by read position. Substitution types are collapsed based on the
        reference or expected base, with only six substitution types being reported:
        `A>C`, `A>G`, `A>T`, `C>A`, `C>G` and `C>T`.
        For example, `T>G` is grouped in with `A>C`.

        The following are reads / bases are excluded from the analysis:

        * Unmapped reads
        * Reads marked as failing vendor quality
        * Reads marked as duplicates (unless `--include-duplicates` was specified)
        * Secondary and supplemental records
        * Soft-clipped bases in records
        * Reads with MAPQ < `--min-mapping-quality` (default: `20`)
        * Bases with base quality < `--min-base-quality` (default: `0`)
        * Bases where either the read base or the reference base is non-ACGT
        ''',
        plot = linegraph.plot(linegraph_data, pconfig)
    )

    # Add to general stats table
    headers = {
        'error_rate': {
            'title': '% Error',
            'description': 'Percent error across all read positions',
            'min': 0,
            'max': 100.0,
            'scale': 'GnYlRd',
            'suffix': '%',
            'format': '{:,.2f}',
            'modify': lambda x: 100.0 * x
        }
    }
    self.general_stats_addcols(error_rates, headers)

    return len(all_data)
