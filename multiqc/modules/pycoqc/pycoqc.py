#!/usr/bin/env python

""" MultiQC module to parse output from pycoQC """

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table, linegraph
from multiqc import config
from collections import OrderedDict
import yaml
import logging

log = logging.getLogger(__name__)




class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='pycoQC', anchor='pycoqc',
        href="https://a-slide.github.io/pycoQC/",
        info="computes metrics and generates interactive QC plots for Oxford Nanopore technologies sequencing data")
        self.pycoqc_data = {}
        for f in self.find_log_files('pycoqc'):
            if f['s_name'] in self.pycoqc_data:
                log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], f['s_name']))
            self.pycoqc_data[f['s_name']] = self.load_data(f['f'])

        if len(self.pycoqc_data) == 0:
            raise UserWarning

        self.table_data, self.reads_data, self.bases_data, self.read_length_plot_data, self.quality_plot_data = self.parse_data()

        self.general_stats_headers = self.setup_stats_header()
        self.general_stats_addcols(self.table_data, self.general_stats_headers)
        self.write_data_file(self.table_data, 'multiqc_pycoqc')

        self.pycoqc_table_headers = self.setup_pycoqc_table_headers()
        self.add_section (
            name = 'Statistics',
            anchor = 'pycoqc_stats',
            plot = table.plot(self.table_data, self.pycoqc_table_headers)
        )

        self.read_bar_config = self.setup_read_bar_config()
        cats = ['passed_reads', 'non_passed_reads']
        self.add_section(
            name = 'Reads counts',
            anchor = 'pycoqc_reads',
            description = 'Number of sequenced reads passing / failing the QC thresholds.',
            plot = bargraph.plot(self.reads_data, cats, self.read_bar_config)
        )

        self.bases_bar_config = self.setup_bases_bar_config()
        cats = ['passed_bases', 'non_passed_bases']
        self.add_section(
            name = 'Bases counts',
            anchor = 'pycoqc_bases',
            description = 'Number of sequenced bases passing / failing the QC thresholds.',
            plot = bargraph.plot(self.bases_data, cats, self.bases_bar_config)
        )

        if self.read_length_plot_data[0]:
            self.read_length_config = self.setup_read_length_config()
            self.add_section(
                name = 'Read length',
                anchor = 'pycoqc_read_len',
                description = 'Distribution of read lengt for all / passed reads.',
                plot = linegraph.plot(self.read_length_plot_data, self.read_length_config)
            )

        if self.quality_plot_data[0]:
            self.qual_config = self.setup_qual_config()
            self.add_section(
                name = 'Quality scores',
                anchor = 'pycoqc_read_qual',
                description = 'Distribution of quality scores for all / passed reads.',
                plot = linegraph.plot(self.quality_plot_data, self.qual_config)
            )

        self.run_duration_headers = self.setup_run_duration_headers()
        self.add_section(
            name = 'Run Duration',
            anchor = 'pycoqc_run_duration',
            description = 'Run Duration',
            plot = table.plot(self.table_data, self.run_duration_headers)
        )

    def load_data(self, f):
        data = {}
        try:
            data = yaml.load(f, Loader=yaml.FullLoader)
        except Exception as e:
            log.warn("An error occurred when opening '{}'.".format(f))
        return data

    def parse_data(self):
        data_for_table = {}
        reads_data = {}
        bases_data = {}
        length_plot_all = {}
        length_plot_pass = {}
        qual_plot_all = {}
        qual_plot_pass = {}

        for sample, sample_data in self.pycoqc_data.items():
            data_for_table[sample] = {
                'all_median_read_length': sample_data['All Reads']['basecall']['len_percentiles'][50],
                'all_median_phred_score': sample_data['All Reads']['basecall']['qual_score_percentiles'][50],
                'all_n50': sample_data['All Reads']['basecall']['N50'],
                'all_run_duration': sample_data['All Reads']['run']['run_duration'],
                'all_channels': sample_data['All Reads']['run']['active_channels'],
                'all_reads': sample_data['All Reads']['basecall']['reads_number'],
                'all_bases': sample_data['All Reads']['basecall']['bases_number'],
                'passed_median_read_length': sample_data['Pass Reads']['basecall']['len_percentiles'][50],
                'passed_median_phred_score': sample_data['Pass Reads']['basecall']['qual_score_percentiles'][50],
                'passed_n50': sample_data['Pass Reads']['basecall']['N50'],
                'passed_channels': sample_data['Pass Reads']['run']['active_channels'],
                'passed_reads': sample_data['Pass Reads']['basecall']['reads_number'],
                'passed_bases': sample_data['Pass Reads']['basecall']['bases_number'],
                }

            reads_data[sample] = {
                'passed_reads': sample_data['Pass Reads']['basecall']['reads_number'],
                'non_passed_reads': sample_data['All Reads']['basecall']['reads_number'] - sample_data['Pass Reads']['basecall']['reads_number'],
            }

            bases_data[sample] = {
                'passed_bases': sample_data['Pass Reads']['basecall']['bases_number'],
                'non_passed_bases': sample_data['All Reads']['basecall']['bases_number'] - sample_data['Pass Reads']['basecall']['bases_number'],
            }
            try:
                length_x_vals_all = sample_data['All Reads']['basecall']['len_hist']['x']
                length_y_vals_all = sample_data['All Reads']['basecall']['len_hist']['y']
                length_plot_all[sample] = dict(zip(length_x_vals_all,length_y_vals_all))
                length_x_vals_pass = sample_data['Pass Reads']['basecall']['len_hist']['x']
                length_y_vals_pass = sample_data['Pass Reads']['basecall']['len_hist']['y']
                length_plot_pass[sample] = dict(zip(length_x_vals_pass,length_y_vals_pass))

                qual_x_vals_all = sample_data['All Reads']['basecall']['qual_score_hist']['x']
                qual_y_vals_all = sample_data['All Reads']['basecall']['qual_score_hist']['y']
                qual_plot_all[sample] = dict(zip(qual_x_vals_all,qual_y_vals_all))
                qual_x_vals_pass = sample_data['Pass Reads']['basecall']['qual_score_hist']['x']
                qual_y_vals_pass = sample_data['Pass Reads']['basecall']['qual_score_hist']['y']
                qual_plot_pass[sample] = dict(zip(qual_x_vals_pass,qual_y_vals_pass))
            except KeyError:
                log.warn("No plot data found for sample '{}'. Please make sure you are using pycoQC v2.5.0.20 or newer.".format(sample))
        return data_for_table, reads_data, bases_data, [length_plot_pass, length_plot_all], [qual_plot_pass, qual_plot_all]

    def setup_stats_header(self):
        general_stats_headers = OrderedDict()
        general_stats_headers['all_median_read_length'] = {
            'title': 'Median Read Length (All)',
            'description': 'Median Read Length all',
            'scale': 'BuPu',
            'shared_key': 'median_read_len',
        }
        general_stats_headers['passed_median_read_length'] = {
            'title': 'Median Read Length (Pass)',
            'description': 'Median Read Length pass',
            'scale': 'BuPu',
            'shared_key': 'median_read_len',
        }
        general_stats_headers['all_reads'] = {
            'title': '{} reads (All)'.format(config.long_read_count_prefix),
            'description': 'Number of reads all ({})'.format(config.long_read_count_desc),
            'scale': 'BuGn',
            'modify': lambda x: x * config.long_read_count_multiplier,
            'shared_key': 'long_read_count',
        }
        general_stats_headers['passed_reads'] = {
            'title': '{} reads (Pass)'.format(config.long_read_count_prefix),
            'description': 'Number of reads pass ({})'.format(config.long_read_count_desc),
            'scale': 'BuGn',
            'modify': lambda x: x * config.long_read_count_multiplier,
            'shared_key': 'long_read_count',
        }
        general_stats_headers['all_bases'] = {
            'title': '{} Bases (All)'.format(config.base_count_prefix),
            'description': 'Number of bases all ({})'.format(config.base_count_desc),
            'scale': 'OrRd',
            'modify': lambda x: x * config.base_count_multiplier,
            'shared_key': 'base_count',
        }
        general_stats_headers['passed_bases'] = {
            'title': '{} Bases (Pass)'.format(config.base_count_prefix),
            'description': 'Number of bases pass ({})'.format(config.base_count_desc),
            'scale': 'OrRd',
            'modify': lambda x: x * config.base_count_multiplier,
            'shared_key': 'base_count',
        }
        return general_stats_headers

    def setup_pycoqc_table_headers(self):
        pycoqc_table_headers = OrderedDict()
        pycoqc_table_headers['all_n50'] = {
            'namespace': 'pycoQC',
            'title': 'N50 (All)',
            'description': 'N50',
            'scale': 'BuGn',
            'shared_key': 'n50',
        }
        pycoqc_table_headers['passed_n50'] = {
            'namespace': 'pycoQC',
            'title': 'N50 (Pass)',
            'description': 'N50',
            'scale': 'BuGn',
            'shared_key': 'n50',
        }
        pycoqc_table_headers['all_median_phred_score'] = {
            'namespace': 'pycoQC',
            'title': 'Median PHRED score (All)',
            'description': 'Median PHRED score',
            'scale': 'PuRd',
            'shared_key': 'phred',
        }
        pycoqc_table_headers['passed_median_phred_score'] = {
            'namespace': 'pycoQC',
            'title': 'Median PHRED score (Pass)',
            'description': 'Median PHRED score',
            'scale': 'PuRd',
            'shared_key': 'phred',
        }
        pycoqc_table_headers['all_channels'] = {
            'namespace': 'pycoQC',
            'title': 'Active Channels (All)',
            'description': 'Number of active channels',
            'scale': 'OrRd',
            'shared_key': 'channels',
        }
        pycoqc_table_headers['passed_channels'] = {
            'namespace': 'pycoQC',
            'title': 'Active Channels (Pass)',
            'description': 'Number of active channels',
            'scale': 'OrRd',
            'shared_key': 'channels',
        }
        return pycoqc_table_headers

    def setup_run_duration_headers(self):
        run_duration_headers = OrderedDict()
        run_duration_headers['all_run_duration'] = {
            'namespace': 'pycoQC',
            'title': 'Run duration',
            'description': 'Run duration (h)',
            'scale': 'PuBuGn'
        }
        return(run_duration_headers)

    def setup_read_bar_config(self):
        read_length_bar_plot_config = {
            'id': 'pycoqc_length_plot',
            'title': 'pycoQC: Read Length',
            'ylab': 'Sample',
        }
        return read_length_bar_plot_config

    def setup_bases_bar_config(self):
        read_qual_bar_plot_config = {
            'id': 'pycoqc_qual_plot',
            'title': 'pycoQC: Read Quality',
            'ylab': 'Sample',
        }
        return read_qual_bar_plot_config

    def setup_read_length_config(self):
        read_length_config = {
            'data_labels': [
                {'name': 'Pass', 'ylab': 'Read Density', 'xlab': 'Basecalled Length'},
                {'name': 'All', 'ylab': 'Read Density', 'xlab': 'Basecalled Length'}
            ]
        }
        return read_length_config

    def setup_qual_config(self):
        qual_config = {
            'data_labels': [
                {'name': 'Pass', 'ylab': 'Read Density', 'xlab': 'Read Quality Score'},
                {'name': 'All', 'ylab': 'Read Density', 'xlab': 'Read Quality Score'}
            ]
        }
        return qual_config
