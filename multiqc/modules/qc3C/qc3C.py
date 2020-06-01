from __future__ import print_function
from multiqc.modules.base_module import BaseMultiqcModule
from _collections import OrderedDict
from multiqc import config
from multiqc.plots import table, bargraph
from os.path import dirname, basename
import numpy as np
import json
import re
import logging

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='qc3C', anchor='qc3C',
                                            href="http://github.com/cerebis/qc3C",
                                            info="provides reference-free and BAM based quality control for Hi-C data")

        self.qc3c_data = {}
        self.digest_junctions = {}

        for f in self.find_log_files('qc3C', filehandles=True):
            self.parse_qc3c_log(f)

        self.qc3c_data = self.ignore_samples(self.qc3c_data)

        if len(self.qc3c_data) == 0:
            raise UserWarning
        log.info("Found {} reports".format(len(self.qc3c_data)))

        self.write_data_file(self.qc3c_data, 'multiqc_qc3c')

        self.add_section(
            name='qc3C: Runtime Parameters',
            anchor='runtime-parameters',
            plot=self.qc3c_runtime_table())

        self.add_section(
            name='qc3C: Estimated Hi-C fraction',
            anchor='hic-fraction',
            plot=self.qc3c_signal_table())

        self.add_section(
            name='qc3C: Breakdown of Parsed Reads',
            anchor='acceptance_plot',
            plot=self.qc3c_acceptance_plot())

        self.add_section(
            name='qc3C: Putative Junction Content',
            anchor='qc3c_signal_plot',
            plot=self.qc3c_signal_plot())

        # self.add_section(
        #     name='qc3C: Possible junction combinations',
        #     anchor='runtime-parameters',
        #     plot=self.qc3c_combo_table())

        self.add_section(
            name='qc3C: Junction frequency breakdown',
            anchor='qc3c_junction_plot',
            plot=self.qc3c_junction_plot())

    @staticmethod
    def _drop_time(s):
        """ remove time """
        m = re.match(r'(.*) .+$', s)
        return s if m is None else m.group(1)

    @staticmethod
    def _drop_name(s):
        """ drop leading program name from version string """
        return s.split()[-1]

    def qc3c_runtime_table(self):

        config = {'id': 'qc3c_runtime_table', 'namespace': 'qc3C', 'col1_header': 'Sample'}
        headers = OrderedDict({
            'run_timestamp': {'title': 'Date',
                              'description': "Analysis time stamp",
                              'modify': MultiqcModule._drop_time},
            'mode': {'title': 'Run Mode',
                     'description': 'Analysis mode used'},
            'kmer_size': {'title': 'k',
                          'description': 'Library k-mer size',
                          'min': 0, 'format': '{:d}', 'scale': False},
            'enzymes': {'title': 'Digest',
                        'description': 'Enzymes used in digest'},
            'seed': {'title': 'Seed',
                     'description': 'Random seed',
                     'format': '{:d}', 'scale': False},
            # TODO this might not exist if users have not requested a limit
            'max_obs': {'title': 'Max obs',
                        'description': 'User specified maximum number of observations',
                        'min': 0, 'format': '{:d}', 'scale': False},
            'max_freq_quantile': {'title': 'Quantile',
                                  'description': 'Quantile cut-off for low-pass k-mer frequency filter',
                                  'min': 0, 'max': 1, 'format': '{:g}'},
            'sample_rate': {'title': 'Sample rate',
                            'description': 'Sub-sampling probability',
                            'min': 0, 'max': 1, 'format': '{:g}'},
            'mean_insert': {'title': 'Insert size',
                            'description': 'User-specified insert size',
                            'min': 0, 'format': '{:.0f,}', 'suffix': 'bp'},
            'max_freq': {'title': 'Max freq',
                         'description': 'Maximum k-mer frequency after quantile filtering',
                         'min': 0, 'format': '{:d,}'},
            'mean_readlen': {'title': 'Read length',
                             'description': 'Observed average read length',
                             'format': '{:,.0f}', 'suffix': 'bp'},
        })
        return table.plot(self.qc3c_data, headers, config)

    def qc3c_signal_table(self):

        config = {'id': 'qc3c_signal_table', 'namespace': 'qc3C', 'col1_header': 'Sample'}
        headers = OrderedDict({
            'run_timestamp': {'title': 'Date',
                              'description': "Analysis time stamp",
                              'modify': MultiqcModule._drop_time},
            'qc3C_version': {'title': 'Version',
                             'description': 'Program version',
                             'modify': MultiqcModule._drop_name},
            'mode': {'title': 'Run Mode',
                     'description': 'Analysis mode used'},
            'unobs_fraction': {'title': 'Unobserved extent',
                               'description': 'Estimated fraction of total fragment extent that was unobservable',
                               'shared_key': 'unobs_mean',
                               'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Reds'},
            'raw_fraction': {'title': 'Raw Hi-C estimate',
                             'description': 'Raw estimate of Hi-C fraction from observable extent',
                             'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Greens'},
            'adj_fraction': {'title': 'Adjusted Hi-C estimate',
                             'description': 'Estimate of Hi-C fraction adjusted for unobserved extent',
                             'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Greens'},
        })
        return table.plot(self.qc3c_data, headers, config)

    # def qc3c_combo_table(self):
    #
    #     config = {'id': 'qc3c_combo_table', 'namespace': 'qc3C', 'col1_header': 'Sample'}
    #     headers = OrderedDict({
    #         'enz5p': {'title': '5-prime end',
    #                   'description': '5-prime side of ligation junction'},
    #         'enz3p': {'title': '5-prime end',
    #                   'description': '5-prime side of ligation junction'},
    #         'raw_fraction': {'title': 'Raw Hi-C estimate',
    #                          'description': 'Raw estimate of Hi-C fraction from observable extent',
    #                          'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Greens'},
    #         'adj_fraction': {'title': 'Adjusted Hi-C estimate',
    #                          'description': 'Estimate of Hi-C fraction adjusted for unobserved extent',
    #                          'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Greens'},
    #     })
    #     return table.plot(self.qc3c_data, headers, config)
    #
    def qc3c_acceptance_plot(self):
        config = {'id': 'qc3c_acceptance_plot',
                  'namespace': 'qc3C',
                  'ylab': 'Number of Reads',
                  'hide_zero_cats': False,
                  'cpswitch_counts_label': 'Number of Reads',}
        categories = OrderedDict({
            'n_skipped': {'name': 'Skipped', 'color': '#b3b3b3',},
            'n_too_short': {'name': 'Too short', 'color': '#ffd92f',},
            'n_no_flank': {'name': 'No flank', 'color': '#66c2a5',},
            'n_ambiguous': {'name': 'Ambiguous', 'color': '#e78ac3',},
            'n_high_cov': {'name': 'High cov', 'color': '#8da0cb',},
            'n_accepted_reads': {'name': 'Accepted', 'color': '#fc8d62',},
        })
        return bargraph.plot(self.qc3c_data, categories, config)

    def qc3c_signal_plot(self):
        config = {'id': 'qc3c_signal_plot',
                  'namespace': 'qc3C',
                  'ylab': 'Number of Reads',
                  'hide_zero_cats': False,
                  'cpswitch_counts_label': 'Number of Reads',}
        categories = OrderedDict({
            'n_without_junc': {'name': 'Without junc', 'color': '#b3b3b3'},
            'n_with_junc': {'name': 'With junc', 'color': '#fc8d62'},
            # 'n_cs_start': {'name': 'Cut-site start'},
        })
        return bargraph.plot(self.qc3c_data, categories, config)

    def qc3c_junction_plot(self):
        config = {'id': 'qc3c_frequency_plot',
                  'namespace': 'qc3C',
                  'hide_zero_cats': False,
                  'use_legend': False,}
        categories = OrderedDict({vi: {'name': vi} for v in self.digest_junctions.values() for vi in v})
        return bargraph.plot(self.qc3c_data, categories, config)

    def parse_qc3c_log(self, f):
        try:
            parsed = json.load(f['f'])
        except json.JSONDecodeError:
            log.warning("Could not parse qc3C JSON: '{}'".format(f['fn']))
            return None

        s_name = self.clean_s_name(basename(f['root']), dirname(f['root']))
        if parsed is not None:
            if s_name in self.qc3c_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))

        for k in 'raw_fraction', 'adj_fraction', 'unobs_fraction':
            parsed[k] = np.array(parsed[k]).mean() * 100

        self.qc3c_data[s_name] = {'qc3C_version': parsed['runtime_info']['qc3C_version'],
                                  'run_timestamp': parsed['runtime_info']['run_timestamp'],
                                  'mode': parsed['mode'],
                                  # TODO this must be conditional -- check all such expectations
                                  'kmer_size': parsed['input_args']['kmer_size'],
                                  'enzymes': ', '.join(parsed['input_args']['enzymes']),
                                  'seed': parsed['input_args']['seed'],
                                  'sample_rate': parsed['input_args']['sample_rate'],
                                  'max_freq': parsed['input_args']['max_coverage'],
                                  'mean_insert': parsed['input_args']['mean_insert'],
                                  'max_freq_quantile': parsed['input_args']['max_freq_quantile'],
                                  'max_obs': parsed['input_args']['max_obs'],
                                  'n_skipped': parsed['n_parsed_reads'] - parsed['n_analysed_reads'],
                                  'n_analysed_reads': parsed['n_analysed_reads'],
                                  'n_too_short': parsed['n_too_short'],
                                  'n_no_flank': parsed['n_no_flank'],
                                  'n_ambiguous': parsed['n_ambiguous'],
                                  'n_high_cov': parsed['n_high_cov'],
                                  'n_accepted_reads': parsed['n_accepted_reads'],
                                  'n_without_junc': parsed['n_without_junc'],
                                  'n_with_junc': parsed['n_with_junc'],
                                  'mean_readlen': parsed['mean_readlen'],
                                  'n_cs_start': parsed['cs_start'],
                                  'raw_fraction': parsed['raw_fraction'],
                                  'adj_fraction': parsed['adj_fraction'],
                                  'unobs_fraction': parsed['unobs_fraction'],
                                  }

        # include the junction frequencies (1 or many depending on digest)
        self.qc3c_data[s_name].update(parsed['junction_frequency'])

        # retain the junction keys for plotting, as these may vary per experiment
        self.digest_junctions[s_name] = sorted(parsed['junction_frequency'].keys())

        self.add_data_source(f, s_name)
