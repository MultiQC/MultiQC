from __future__ import print_function
from multiqc.modules.base_module import BaseMultiqcModule
from collections import OrderedDict, defaultdict
from multiqc import config
from multiqc.plots import table, bargraph, linegraph
from os.path import dirname, basename
import numpy as np
import json
import re
import logging

log = logging.getLogger(__name__)


#
# A 64 color palette defined in four gradient blocks. Red->Yellow, Yellow->Green, Green->Blue, Blue->Deep Blue
# The method color_picker uses these to sample a gradient with varying sparseness

# red to yellow
ry_1 = ["#EF476F", "#F0506E", "#F1586E", "#F2616D", "#F36A6D", "#F4726C", "#F57B6C", "#F6836B", "#F78C6B",
        "#F8956A", "#F99D69", "#FAA669", "#FBAF68", "#FCB768", "#FDC067", "#FEC867", "#FFD166"]

# yellow to green
yg_2 = ["#FFD166", "#EFD16A", "#E0D26D", "#D0D271", "#C1D275", "#B1D378", "#A2D37C", "#92D37F", "#83D483",
        "#73D487", "#63D48A", "#54D48E", "#44D592", "#35D595", "#25D599", "#16D69C", "#06D6A0"]

# green to blue
gb_3 = ["#06D6A0", "#07D1A1", "#07CDA2", "#08C8A3", "#09C3A5", "#09BEA6", "#0ABAA7", "#0BB5A8", "#0CB0A9",
        "#0CABAA", "#0DA7AB", "#0EA2AC", "#0E9DAE", "#0F98AF", "#1094B0", "#108FB1", "#118AB2"]

# blue to dark blue
bd_4 = ["#118AB2", "#1085AB", "#107FA4", "#0F7A9E", "#0E7597", "#0E7090", "#0D6A89", "#0C6582", "#0C607C",
        "#0B5B75", "#0A556E", "#0A5067", "#094B60", "#08465A", "#084053", "#073B4C"]

# the 64-color gradient
grad64 = ry_1[:-1] + yg_2[:-1] + gb_3[:-1] + bd_4


def color_picker(degen):
    """
    Select colours from a 64 colour gradient, where an effort is made to use the full
    width of the gradient depending on the number of elements required.
    :param degen: a list of degeneracies (number of colours per object)
    :return: a flat list of colours equal the sum of degen
    """
    if len(degen) == 1:
        # a single non-ambiguous enzyme, lets make this blue
        if degen[0] not in {1, 4, 16}:
            raise ValueError(f'got {degen[0]} junc_degen values can only be 1, 4 or 16')
        return grad64[0::64 // degen[0]]
    else:
        cols = []
        for n, jd in enumerate(degen):
            if jd not in {1, 4, 16}:
                raise ValueError(f'got {jd} when junc_degen values can only be 1, 4 or 16')
            cols += grad64[16*n:16*n+16:16//jd]
        return cols


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='qc3C', anchor='qc3C',
                                            href="http://github.com/cerebis/qc3C",
                                            info="provides reference-free and BAM based quality control for Hi-C data")

        self.qc3c_data = defaultdict(dict)
        self.digest_junctions = defaultdict(dict)

        for f in self.find_log_files('qc3C', filehandles=True):
            self.parse_qc3c_log(f)

        for k in self.qc3c_data:
            self.qc3c_data[k] = self.ignore_samples(self.qc3c_data[k])

        # check that any non-empty records were stored under one of the two analysis modes
        n_reports = len(self.qc3c_data['kmer']) + len(self.qc3c_data['bam'])
        if n_reports == 0:
            raise UserWarning('No reports found')
        log.info("Found {} reports".format(n_reports))

        if len(self.qc3c_data['bam']) > 0:

            self.write_data_file(self.qc3c_data['bam'], 'multiqc_qc3c_bam')
            log.info('Found {} BAM analysis reports'.format(len(self.qc3c_data['bam'])))

            self.add_section(
                name='qc3C BAM: Runtime details',
                anchor='bam-runtime-parameters',
                plot=self.bam_runtime_table())

            self.add_section(
                name='qc3C BAM: Pair breakdown',
                anchor='bam-hic-fraction',
                plot=self.bam_signal_table())

            self.add_section(
                name='qc3C BAM: Long-range pairs',
                anchor='bam-longrange-plot',
                plot=self.bam_longrange_plot())

            self.add_section(
                name='qc3C BAM: Read parsing',
                anchor='bam-acceptance-plot',
                plot=self.bam_acceptance_plot())

            self.add_section(
                name='qc3C BAM: Pair validation',
                anchor='bam-valid-plot',
                plot=self.bam_valid_plot())

            self.add_section(
                name='qc3C BAM: Expected junctions',
                anchor='bam-junction-plot',
                plot=self.bam_junction_plot())

            self.add_section(
                name='qc3C BAM: Distribution of fragment separation',
                anchor='bam-fragment-histogram',
                plot=self.bam_fragment_histogram())

        if len(self.qc3c_data['kmer']) > 0:

            self.write_data_file(self.qc3c_data['kmer'], 'multiqc_qc3c_kmer')
            log.info('Found {} k-mer analysis reports'.format(len(self.qc3c_data['kmer'])))

            self.add_section(
                name='qc3C K-mer: Runtime details',
                anchor='kmer-runtime-parameters',
                plot=self.kmer_runtime_table())

            self.add_section(
                name='qc3C K-mer: Estimated Hi-C fraction',
                anchor='kmer-hic-fraction',
                plot=self.kmer_signal_table())

            self.add_section(
                name='qc3C K-mer: Breakdown of parsed reads',
                anchor='kmer-acceptance-plot',
                plot=self.kmer_acceptance_plot())

            self.add_section(
                name='qc3C K-mer: Putative junction content',
                anchor='kmer-signal-plot',
                plot=self.kmer_signal_plot())

            self.add_section(
                name='qc3C K-mer: Junction frequency breakdown',
                anchor='kmer-junction-plot',
                plot=self.kmer_junction_plot())

    @staticmethod
    def _drop_time(s):
        """ remove time """
        m = re.match(r'(.*) .+$', s)
        return s if m is None else m.group(1)

    @staticmethod
    def _drop_name(s):
        """ drop leading program name from version string """
        return s.split()[-1]

    def bam_runtime_table(self):

        config = {'id': 'bam_runtime_table', 'namespace': 'qc3C', 'col1_header': 'Sample'}
        headers = OrderedDict({
            'run_timestamp': {'title': 'Date',
                              'description': "Analysis time stamp",
                              'modify': MultiqcModule._drop_time},
            'mode': {'title': 'Run Mode',
                     'description': 'Analysis mode used'},
            'min_mapq': {'title': 'Min MapQ',
                         'description': 'Minimum accepted mapping quality',
                        'min': 0, 'format': '{:d}', 'scale': False},
            'enzymes': {'title': 'Digest',
                        'description': 'Enzymes used in digest'},
            'seed': {'title': 'Seed',
                     'description': 'Random seed',
                     'format': '{:d}', 'scale': False},
            'max_obs': {'title': 'Max obs',
                        'description': 'User specified maximum number of observations',
                        'min': 0, 'format': '{:,d}', 'scale': False, 'modify': lambda x: 'n/a' if x == -1 else x},
            'n_accepted_pairs': {'title': 'Accepted pairs',
                                 'description': 'Number of pairs accepted for analysis',
                                 'min': 0, 'format': '{:,d}'},
            'sample_rate': {'title': 'Sample rate',
                            'description': 'Sub-sampling probability',
                            'min': 0, 'max': 1, 'format': '{:g}'},
            'obs_insert_mean': {'title': 'Insert mean',
                                'description': 'Estimated mean insert size',
                                'min': 0, 'format': '{:,.0f}', 'suffix': 'bp'},
            'obs_insert_median': {'title': 'Insert median',
                                  'description': 'Estimated median insert size',
                                  'min': 0, 'format': '{:,.0f}', 'suffix': 'bp'},
            'mean_readlen': {'title': 'Read length',
                             'description': 'Observed average read length',
                             'format': '{:,.0f}', 'suffix': 'bp'},
        })
        return table.plot(self.qc3c_data['bam'], headers, config)

    def bam_longrange_plot(self):
        config = {'id': 'bam_longrange_plot',
                  'namespace': 'qc3C',
                  'ylab': 'Number of Reads',
                  'hide_zero_cats': False,
                  'cpswitch_counts_label': 'Number of Reads',}

        categories = OrderedDict({
            'n_1kb_pairs': {'name': '>1000 bp',},
            'n_5kb_pairs': {'name': '>5000 bp',},
            'n_10kb_pairs': {'name': '>10000 bp',},
            'n_accepted_pairs': {'name': 'Accepted',}
        })

        return bargraph.plot(self.qc3c_data['bam'], categories, config)

    def bam_acceptance_plot(self):
        config = {'id': 'bam_acceptance_plot',
                  'namespace': 'qc3C',
                  'ylab': 'Number of Reads',
                  'hide_zero_cats': False,
                  'cpswitch_counts_label': 'Number of Reads',}
        categories = OrderedDict({
            'n_skipped_reads': {'name': 'Skipped',},
            'n_unmapped_reads': {'name': 'Unmapped',},
            'n_low_mapq_reads': {'name': 'Low mapq',},
            'n_secondary_reads': {'name': 'Secondary',},
            'n_supplementary_reads': {'name': 'Supplementary',},
            'n_weak_mapping_reads': {'name': 'Weak mapping',},
            'n_ref_term_reads': {'name': 'Truncated',},
            'n_accepted_reads': {'name': 'Accepted',},
        })
        return bargraph.plot(self.qc3c_data['bam'], categories, config)

    def bam_signal_table(self):

        config = {'id': 'bam_signal_table',
                  'namespace': 'qc3C',
                  'hide_zero_cats': False,
                  'col1_header': 'Sample'}
        headers = OrderedDict({
            'run_timestamp': {'title': 'Date',
                              'description': "Analysis time stamp",
                              'modify': MultiqcModule._drop_time},
            'mode': {'title': 'Run Mode',
                     'description': 'Analysis mode used'},
            'p_trans_pairs': {'title': 'Trans pairs',
                              'description': 'Fraction of inter-contig pairs',
                              'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Greens'},
            'p_cis_pairs': {'title': 'Cis pairs',
                            'description': 'Fraction of intra-contig pairs',
                            'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Greens'},
            'p_fully_aligned': {'title': 'Fully aligned',
                                'description': 'Fraction of full alignments',
                                'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Greens'},
            'p_align_term': {'title': 'Trunc aligned',
                             'description': 'Fraction of truncated alignments',
                             'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Greens'},
            'p_no_site_end':{'title': 'No site end',
                             'description': 'Fraction alignments not ending in a cutsite',
                             'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Greens'},
            'p_short_inserts':{'title': 'Short insert',
                               'description': 'Fraction small-separation pairs < 1000bp',
                               'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Reds'},
            'unobs_fraction': {'title': 'Unobserved extent',
                               'description': 'Estimated fraction of total fragment extent that was unobservable',
                               'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Reds'},
            'p_informative_fr': {'title': "Valid FR",
                                 'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Greens'},
            'p_informative_rf': {'title': "Valid RF",
                                 'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Greens'},
            'p_informative_ffrr': {'title': "Valid FF|RR",
                                    'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Greens'},
            'p_uninformative_religation': {'title': "Religation",
                                           'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Reds'},
            'p_uninformative_dangling_ends': {'title': "Dangling End",
                                              'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Reds'},
            'p_uninformative_self_circle': {'title': "Self-circle",
                                            'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Reds'},
            'p_uninformative_ffrr': {'title': "Invalid FF|RR",
                                     'min': 0, 'max': 100, 'suffix': '%', 'scale': 'Reds'},

        })
        return table.plot(self.qc3c_data['bam'], headers, config)

    def bam_valid_plot(self):
        config = {'id': 'bam_valid_plot',
                  'namespace': 'qc3C',
                  'ylab': 'Number of Reads',
                  'hide_zero_cats': False,
                  'cpswitch_counts_label': 'Number of Reads',}
        categories = OrderedDict({
            'n_informative_fr': {'name': "Valid FR", 'color': '#a1d99b'},
            'n_informative_rf': {'name': "Valid RF", 'color': '#74c476'},
            'n_informative_ffrr': {'name': "Valid FF|RR", 'color': '#41ab5d'},
            'n_uninformative_religation': {'name': "Religation", 'color': '#fcbba1'},
            'n_uninformative_dangling_ends': {'name': "Dangling End", 'color': '#fc9272'},
            'n_uninformative_self_circle': {'name': "Self-circle", 'color': '#fb6a4a'},
            'n_uninformative_ffrr': {'name': "Invalid FF|RR", 'color': '#ef3b2c'},
        })
        return bargraph.plot(self.qc3c_data['bam'], categories, config)

    def bam_junction_plot(self):
        config = {'id': 'bam_junction_plot',
                  'namespace': 'qc3C',
                  'hide_zero_cats': False,
                  'use_legend': False,}

        categories = OrderedDict()
        for v in self.digest_junctions['bam'].values():
            for vi in v:
                categories[vi['name']] = vi

        return bargraph.plot(self.qc3c_data['bam'], categories, config)

    def bam_fragment_histogram(self):

        median_lines = []
        for smpl in self.qc3c_data['bam']:
            median_lines.append({'value': self.qc3c_data['bam'][smpl]['obs_insert_median'],
                                 'color': '#D8E2DC',
                                 'width': 2,
                                 'dashStyle': 'ShortDashDot'
                                 })

        config = {'id': 'bam_fragment_histogram',
                  'namespace': 'qc3C',
                  'xLog': True,
                  'xPlotLines': median_lines,
                  'xlab': 'log10 [Separation] (bp)',
                  'ylab': 'Density',
                  'tt_label': 'x:{point.x:.1f}, y:{point.y:.4f}'
                  }

        data = {}
        for smpl in self.qc3c_data['bam']:
            data[smpl] = self.qc3c_data['bam'][smpl]['frag_hist']

        return linegraph.plot(data, config)

    def kmer_runtime_table(self):

        config = {'id': 'kmer_runtime_table', 'namespace': 'qc3C', 'col1_header': 'Sample'}
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
            'max_obs': {'title': 'Max obs',
                        'description': 'User specified maximum number of observations',
                        'min': 0, 'format': '{:,d}', 'scale': False, 'modify': lambda x: 'n/a' if x == -1 else x},
            'n_accepted_reads': {'title': 'Accepted reads',
                                 'description': 'Number of reads accepted for analysis',
                                 'min': 0, 'format': '{:,d}'},
            'max_freq_quantile': {'title': 'Quantile',
                                  'description': 'Quantile cut-off for low-pass k-mer frequency filter',
                                  'min': 0, 'max': 1, 'format': '{:g}'},
            'sample_rate': {'title': 'Sample rate',
                            'description': 'Sub-sampling probability',
                            'min': 0, 'max': 1, 'format': '{:g}'},
            'mean_insert': {'title': 'Insert size',
                            'description': 'User-specified insert size',
                            'min': 0, 'format': '{:,.0f}', 'suffix': 'bp'},
            'max_freq': {'title': 'Max freq',
                         'description': 'Maximum k-mer frequency after quantile filtering',
                         'min': 0, 'format': '{:,d}'},
            'mean_readlen': {'title': 'Read length',
                             'description': 'Observed average read length',
                             'format': '{:,.0f}', 'suffix': 'bp'},
        })
        return table.plot(self.qc3c_data['kmer'], headers, config)

    def kmer_signal_table(self):

        config = {'id': 'kmer_signal_table', 'namespace': 'qc3C', 'col1_header': 'Sample'}
        headers = OrderedDict({
            'run_timestamp': {'title': 'Date',
                              'description': "Analysis time stamp",
                              'modify': MultiqcModule._drop_time},
            'mode': {'title': 'Run Mode',
                     'description': 'Analysis mode used'},
            'kmer_size': {'title': 'k',
                          'description': 'Library k-mer size',
                          'min': 0, 'format': '{:d}', 'scale': False},
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
        return table.plot(self.qc3c_data['kmer'], headers, config)

    def kmer_acceptance_plot(self):
        config = {'id': 'kmer_acceptance_plot',
                  'namespace': 'qc3C',
                  'ylab': 'Number of Reads',
                  'hide_zero_cats': False,
                  'cpswitch_counts_label': 'Number of Reads',}
        categories = OrderedDict({
            'n_skipped': {'name': 'Skipped', 'color': '#fcbba1'},
            'n_too_short': {'name': 'Too short', 'color': '#fc9272'},
            'n_no_flank': {'name': 'No flank', 'color': '#fb6a4a'},
            'n_ambiguous': {'name': 'Ambiguous', 'color': '#ef3b2c'},
            'n_high_cov': {'name': 'High cov', 'color': '#cb181d'},
            'n_accepted_reads': {'name': 'Accepted', 'color': '#41ab5d'},
        })
        return bargraph.plot(self.qc3c_data['kmer'], categories, config)

    def kmer_signal_plot(self):
        config = {'id': 'kmer_signal_plot',
                  'namespace': 'qc3C',
                  'ylab': 'Number of Reads',
                  'hide_zero_cats': False,
                  'cpswitch_counts_label': 'Number of Reads',}
        categories = OrderedDict({
            'n_without_junc': {'name': 'Without junc', 'color': '#ef3b2c'},
            'n_with_junc': {'name': 'With junc', 'color': '#41ab5d'},
        })
        return bargraph.plot(self.qc3c_data['kmer'], categories, config)

    def kmer_junction_plot(self):
        config = {'id': 'kmer_frequency_plot',
                  'namespace': 'qc3C',
                  'hide_zero_cats': False,
                  'use_legend': False,}

        categories = OrderedDict()
        for _cat in self.digest_junctions['kmer'].values():
            for vi in _cat:
                categories[vi['name']] = vi

        return bargraph.plot(self.qc3c_data['kmer'], categories, config)

    def parse_qc3c_log(self, f):

        def _none_to(x, y):
            return y if x is None else y

        try:
            parsed = json.load(f['f'])
        except json.JSONDecodeError:
            log.warning("Could not parse qc3C JSON: '{}'".format(f['fn']))
            return None

        s_name = self.clean_s_name(basename(f['root']), dirname(f['root']))
        if parsed is not None:
            if s_name in self.qc3c_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))

        analysis_mode = parsed['mode']

        if analysis_mode == 'bam':

            # set some variables to shorten the lines below
            inf = parsed['classification']['informative']
            uninf = parsed['classification']['uninformative']
            n_cis_pairs = parsed['n_cis_pairs']

            self.qc3c_data['bam'][s_name] = {'qc3C_version': parsed['runtime_info']['qc3C_version'],
                                             'run_timestamp': parsed['runtime_info']['run_timestamp'],
                                             'mode': parsed['mode'],
                                             'enzymes': ', '.join(parsed['input_args']['enzymes']),
                                             'seed': parsed['input_args']['seed'],
                                             'sample_rate': _none_to(parsed['input_args']['sample_rate'], 1),
                                             'max_obs': _none_to(parsed['input_args']['max_obs'], -1),
                                             'n_skipped_reads': parsed['n_skipped_reads'],
                                             'n_unmapped_reads': parsed['n_unmapped'],
                                             'n_analysed_reads': parsed['n_analysed_reads'],
                                             'n_low_mapq_reads': parsed['n_low_mapq'],
                                             'n_ref_len_reads': parsed['n_ref_len'],
                                             'n_secondary_reads': parsed['n_secondary'],
                                             'n_supplementary_reads': parsed['n_supplementary'],
                                             'n_weak_mapping_reads': parsed['n_weak_mapping'],
                                             'n_ref_term_reads': parsed['n_ref_term'],
                                             'n_accepted_reads': parsed['n_accepted_reads'],
                                             'obs_insert_mean': parsed['obs_insert_mean'],
                                             'obs_insert_median': parsed['obs_insert_median'],
                                             'mean_readlen': parsed['mean_readlen'],
                                             'n_analysed_pairs': parsed['n_analysed_pairs'],
                                             'n_accepted_pairs': parsed['n_accepted_pairs'],
                                             'p_trans_pairs': parsed['n_trans_pairs'] / parsed['n_accepted_pairs'] * 100,
                                             'p_cis_pairs': n_cis_pairs / parsed['n_accepted_pairs'] * 100,
                                             'p_fully_aligned': parsed['n_fully_aligned'] / parsed['n_accepted_reads'] * 100,
                                             'p_align_term': parsed['n_align_term'] / parsed['n_accepted_reads'] * 100,
                                             'p_no_site_end': parsed['n_no_site_end'] / parsed['n_accepted_reads'] * 100,
                                             'p_short_inserts': parsed['n_short_inserts'] / parsed['n_accepted_pairs'] * 100,
                                             'n_informative_fr': inf['fr'],
                                             'n_informative_rf': inf['rf'],
                                             'n_informative_ffrr': inf['ffrr'],
                                             'n_uninformative_religation': uninf['religation'],
                                             'n_uninformative_dangling_ends': uninf['dangling_ends'],
                                             'n_uninformative_self_circle': uninf['self_circle'],
                                             'n_uninformative_ffrr': uninf['ffrr'],
                                             'p_informative_fr': inf['fr'] / n_cis_pairs * 100,
                                             'p_informative_rf': inf['rf'] / n_cis_pairs * 100,
                                             'p_informative_ffrr': inf['ffrr'] / n_cis_pairs * 100,
                                             'p_uninformative_religation': uninf['religation'] / n_cis_pairs * 100,
                                             'p_uninformative_dangling_ends': uninf['dangling_ends'] / n_cis_pairs * 100,
                                             'p_uninformative_self_circle': uninf['self_circle'] / n_cis_pairs * 100,
                                             'p_uninformative_ffrr': uninf['ffrr'] / n_cis_pairs * 100,
                                             'n_1kb_pairs': parsed['separation_bins']['counts'][0],
                                             'n_5kb_pairs': parsed['separation_bins']['counts'][1],
                                             'n_10kb_pairs': parsed['separation_bins']['counts'][2],
                                             }

            from itertools import zip_longest
            fhist = {}
            for x, y in zip_longest(parsed['separation_histogram']['mid_points'],
                                    parsed['separation_histogram']['counts']):
                fhist[float(x)] = float(y)
            self.qc3c_data['bam'][s_name]['frag_hist'] = fhist

        elif analysis_mode == 'kmer':

            for k in 'raw_fraction', 'adj_fraction', 'unobs_fraction':
                parsed[k] = np.array(parsed[k]).mean() * 100

            self.qc3c_data['kmer'][s_name] = {'qc3C_version': parsed['runtime_info']['qc3C_version'],
                                              'run_timestamp': parsed['runtime_info']['run_timestamp'],
                                              'mode': parsed['mode'],
                                              'kmer_size': parsed['input_args']['kmer_size'],
                                              'enzymes': ', '.join(parsed['input_args']['enzymes']),
                                              'seed': parsed['input_args']['seed'],
                                              'sample_rate': _none_to(parsed['input_args']['sample_rate'], 1),
                                              'max_freq': parsed['input_args']['max_coverage'],
                                              'mean_insert': parsed['input_args']['mean_insert'],
                                              'max_freq_quantile': parsed['input_args']['max_freq_quantile'],
                                              'max_obs': _none_to(parsed['input_args']['max_obs'], -1),
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
        self.qc3c_data[analysis_mode][s_name].update(parsed['junction_frequency'])

        # calculate the degeneracy of junction sequences per enzymatic combination (5p end =/= 3p end)
        # this can vary due to ambiguous bases in restriction site
        degen_count = {'{}/{}'.format(v['enz5p'], v['enz3p']): 4**v['junction'].count('N')
                       for k, v in parsed['digestion']['junctions'].items()}
        # sort these by "enzyme combo + junction"
        degen_count = OrderedDict(sorted(degen_count.items(), key=lambda x: x[0]))
        # get a palette for this series
        cols = color_picker(list(degen_count.values()))
        # sort these in accordance with that above
        juncs = np.sort(np.array([(k.split(' ')[0], k) for k in parsed['junction_frequency']],
                                 dtype=np.dtype([('a', 'S100'), ('b', 'S100')])))

        # keep a record of how these should be colored per sample
        self.digest_junctions[analysis_mode][s_name] = \
            [{'name': juncs[i][1].decode(), 'color': cols[i]} for i in range(len(juncs))]

        self.add_data_source(f, s_name, section=analysis_mode)
