#!/usr/bin/env python

""" MultiQC module to parse output from Preseq """

from __future__ import print_function
import logging
import numpy as np

from multiqc import config
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Preseq', anchor='preseq',
        href='http://smithlabresearch.org/software/preseq/',
        info="""estimates the complexity of a library, showing how many additional
         unique reads are sequenced for increasing total read count.
         A shallow curve indicates complexity saturation. The dashed line
         shows a perfectly complex library where total reads = unique reads.""")

        # Find and load any Preseq reports
        self.preseq_data = dict()
        self.axis_label = ''
        self.counts_in_1x = None
        for f in self.find_log_files('preseq'):
            parsed_data = self.parse_preseq_logs(f)
            if parsed_data is not None:
                if f['s_name'] in self.preseq_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f)
                self.preseq_data[f['s_name']] = parsed_data

        # Filter to strip out ignored sample names
        self.preseq_data = self.ignore_samples(self.preseq_data)

        if len(self.preseq_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.preseq_data)))

        # Preseq plot
        self.preseq_length_trimmed_plot()


    def parse_preseq_logs(self, f):
        """ Go through log file looking for preseq output """

        lines = f['f'].splitlines()
        header = lines.pop(0)
        using_bases = False
        if header.startswith('TOTAL_READS	EXPECTED_DISTINCT'):
            self.axis_label = 'Molecules'
        elif header.startswith('TOTAL_BASES	EXPECTED_DISTINCT'):
            self.axis_label = 'Bases'
            using_bases = True
        elif header.startswith('total_reads	distinct_reads'):
            self.axis_label = 'Molecules'
        else:
            log.debug("First line of preseq file {} did not look right".format(f['fn']))
            return None

        data = {}
        for l in lines:
            s = l.split()
            # Sometimes the Expected_distinct count drops to 0, not helpful
            if float(s[1]) == 0 and float(s[0]) > 0:
                continue
            data[float(s[0])] = float(s[1])

        # Convert counts to coverage
        read_length = float(getattr(config, 'preseq', {}).get('read_length', 0))
        genome_size = getattr(config, 'preseq', {}).get('genome_size')
        if genome_size is not None:
            try:
                genome_size = float(genome_size)
            except ValueError:
                presets = {'hg19_genome': 2897310462,
                           'hg38_genome': 3049315783,
                           'mm10_genome': 2652783500}
                if genome_size in presets:
                    genome_size = presets[genome_size]
                else:
                    log.warn('The size for genome ' + genome_size + ' is unknown to MultiQC, please specify it '
                             'explicitly or choose one of the following: ' + ', '.join(presets.keys()) +
                             '. Falling back to molecule counts.')
                    genome_size = None
        if genome_size:
            if using_bases:
                self.counts_in_1x = genome_size
            elif read_length:
                self.counts_in_1x = genome_size / read_length
        if self.counts_in_1x:
            data = {k / self.counts_in_1x: v / self.counts_in_1x for k, v in data.items()}
            self.axis_label = 'Coverage'
        return data


    def _read_real_counts(self):
        real_counts_file_raw = None
        real_counts_file_name = None
        for f in self.find_log_files('preseq/real_counts'):
            if real_counts_file_raw is not None:
                log.warn("Multiple Preseq real counts files found, now using {}".format(f['fn']))
            real_counts_file_raw = f['f']
            real_counts_file_name = f['fn']

        real_counts_total = {}
        real_counts_unique = {}
        if real_counts_file_raw is not None:
            try:
                for line in real_counts_file_raw.splitlines():
                    if not line.startswith('#'):
                        cols = line.strip().split()  # Split on any whitespace
                        sn = self.clean_s_name(cols[0], None)
                        if sn in self.preseq_data:
                            if len(cols) >= 2:
                                if cols[1].isdigit():
                                    real_counts_total[sn] = int(cols[1])
                            if len(cols) >= 3:
                                if cols[2].isdigit():
                                    real_counts_unique[sn] = int(cols[2])
            except IOError as e:
                log.error("Error loading real counts file {}: {}".format(real_counts_file_name, str(e)))
            else:
                log.debug("Found {} matching sets of counts from {}".format(len(real_counts_total), real_counts_file_name))

        # Convert real counts to coverage
        if self.counts_in_1x is not None:
            for f, c in real_counts_total.items():
                real_counts_total[f] = float(c) / self.counts_in_1x
            for f, c in real_counts_unique.items():
                real_counts_unique[f] = float(c) / self.counts_in_1x
        return real_counts_total, real_counts_unique


    def _real_counts_to_plot_series(self, real_counts_unique, real_counts_total):
        series = []
        if real_counts_total:
            # Same defaults as HighCharts for consistency
            default_colors = ['#7cb5ec', '#434348', '#90ed7d', '#f7a35c', '#8085e9',
                              '#f15c80', '#e4d354', '#2b908f', '#f45b5b', '#91e8e1']
            for si, sn in enumerate(sorted(self.preseq_data.keys())):
                if sn in real_counts_total:
                    t_reads = float(real_counts_total[sn])
                    point = {
                        'color': default_colors[si % len(default_colors)],
                        'showInLegend': False,
                        'marker': {
                            'enabled': True,
                            'symbol': 'diamond',
                            'lineColor': 'black',
                            'lineWidth': 1,
                        },
                    }
                    if sn in real_counts_unique:
                        u_reads = float(real_counts_unique[sn])
                        point['data'] = [[t_reads, u_reads]]
                        point['name'] = sn + ': actual read count vs. deduplicated read count (externally calculated)'
                        series.append(point)
                        log.debug("Found real counts for {} - Total: {}, Unique: {}"
                                  .format(sn, t_reads, u_reads))
                    else:
                        xvalues = sorted(self.preseq_data[sn].keys())
                        yvalues = sorted(self.preseq_data[sn].values())
                        if t_reads > max(xvalues):
                            log.warning("Total reads for {} ({}) > max preseq value ({}) - "
                                        "skipping this point..".format(sn, t_reads, max(xvalues)))
                        else:
                            interp = np.interp(t_reads, xvalues, yvalues)
                            point['data'] = [[t_reads, interp]]
                            point['name'] = sn + ': actual read count (externally calculated)'
                            series.append(point)
                            log.debug("Found real count for {} - Total: {:.2f} (preseq unique reads: {:.2f})"
                                      .format(sn, t_reads, interp))
        return series


    def preseq_length_trimmed_plot (self):
        """ Generate the preseq plot """

        max_y, sn = max((max(d.values()), s) for s, d in self.preseq_data.items())
        description = ''

        # Plot config
        if self.counts_in_1x is not None:
            precision = '2'
            if max_y > 30:   # no need to be so precise when the depth numbers are high
                precision = '1'
            if max_y > 300:  # when the depth are very high, decimal digits are excessive
                precision = '0'
            tt_label = '<b>{point.x:,.' + precision + 'f}x total</b>: {point.y:,.' + precision + 'f}x unique'
            label_format = '{value}x'
        else:
            tt_label = '<b>{point.x:,.0f} total molecules</b>: {point.y:,.0f} unique molecules'
            label_format = None

        pconfig = {
            'id': 'preseq_plot',
            'title': 'Preseq: Complexity curve',
            'ylab': 'Unique {}'.format(self.axis_label),
            'xlab': 'Total {} (including duplicates)'.format(self.axis_label),
            'ymin': 0,
            'xmin': 0,
            'tt_label': tt_label,
            'yLabelFormat': label_format,
            'xLabelFormat': label_format,
            'extra_series': []
        }

        # Plot the real counts if we have them
        real_counts_total, real_counts_unique = self._read_real_counts()
        if real_counts_unique:
            max_y = max(max_y, max(real_counts_unique.values()))
        pconfig['extra_series'].extend(self._real_counts_to_plot_series(real_counts_unique, real_counts_total))
        if real_counts_unique:
            description += '<p>Points show read count versus deduplicated read counts (externally calculated).</p>'
        elif real_counts_total:
            description += '<p>Points show externally calculated read counts on the curves.</p>'

        # Trim the data to not have a ridiculous x-axis (10Gbp anyone?)
        if getattr(config, 'preseq', {}).get('notrim', False) is not True:
            max_y *= 0.8
            max_x = 0
            for x in sorted(list(self.preseq_data[sn].keys())):
                max_x = max(max_x, x)
                if self.preseq_data[sn][x] > max_y and \
                        x > real_counts_total.get(sn, 0) and x > real_counts_unique.get(sn, 0):
                    break
            pconfig['xmax'] = max_x
            description += "<p>Note that the x axis is trimmed at the point where all the datasets \
                show 80% of their maximum y-value, to avoid ridiculous scales.</p>"

        # Plot perfect library as dashed line
        pconfig['extra_series'].append({
            'name': 'x = y (a perfect library where each read is unique)',
            'data': [[0, 0], [max_y, max_y]],
            'dashStyle': 'Dash',
            'lineWidth': 1,
            'color': '#000000',
            'marker': { 'enabled': False },
            'showInLegend': False,
        })

        self.add_section(
            name = 'Complexity curve',
            description = description,
            anchor = 'preseq-plot',
            plot = linegraph.plot(self.preseq_data, pconfig)
        )
