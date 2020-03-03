#!/usr/bin/env python

""" MultiQC module to parse stats output from pairtools """

import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule


from itertools import combinations, combinations_with_replacement, zip_longest
import os
import re
import json
import numpy as np
from copy import copy, deepcopy
from multiqc.plots import bargraph, linegraph, heatmap
from multiqc.utils import report
from multiqc import config

# Initialise the logger
log = logging.getLogger(__name__)

class ParseError(Exception):
    pass

_SEP = '\t'
_KEY_SEP = '/'

def _from_file(file_handle):
    """create instance of PairCounter from file
    Parameters
    ----------
    file_handle: file handle
    Returns
    -------
    PairCounter
        new PairCounter filled with the contents of the input file
    """
    # fill in from file - file_handle:
    stat_from_file = OrderedDict()
    #
    min_log10_dist=0
    max_log10_dist=9
    log10_dist_bin_step=0.25
    # some variables used for initialization:
    # genomic distance bining for the ++/--/-+/+- distribution
    _dist_bins = (np.r_[0,
        np.round(10**np.arange(min_log10_dist, max_log10_dist+0.001,
                               log10_dist_bin_step))
        .astype(np.int)]
    )

    # establish structure of an empty _stat:
    stat_from_file['total'] = 0
    stat_from_file['total_unmapped'] = 0
    stat_from_file['total_single_sided_mapped'] = 0
    # total_mapped = total_dups + total_nodups
    stat_from_file['total_mapped'] = 0
    stat_from_file['total_dups'] = 0
    stat_from_file['total_nodups'] = 0
    ########################################
    # the rest of stats are based on nodups:
    ########################################
    stat_from_file['cis'] = 0
    stat_from_file['trans'] = 0
    stat_from_file['pair_types'] = {}
    # to be removed:
    stat_from_file['dedup'] = {}

    stat_from_file['cis_1kb+'] = 0
    stat_from_file['cis_2kb+'] = 0
    stat_from_file['cis_4kb+'] = 0
    stat_from_file['cis_10kb+'] = 0
    stat_from_file['cis_20kb+'] = 0
    stat_from_file['cis_40kb+'] = 0

    stat_from_file['chrom_freq'] = OrderedDict()

    stat_from_file['dist_freq'] = OrderedDict([
        ('+-', np.zeros(len(_dist_bins), dtype=np.int)),
        ('-+', np.zeros(len(_dist_bins), dtype=np.int)),
        ('--', np.zeros(len(_dist_bins), dtype=np.int)),
        ('++', np.zeros(len(_dist_bins), dtype=np.int)),
        ])

    # pack distance bins along with dist_freq at least for now
    stat_from_file['dist_bins'] = _dist_bins

    #
    for l in file_handle:
        fields = l.strip().split(_SEP)
        if len(fields) == 0:
            # skip empty lines:
            continue
        if len(fields) != 2:
            # expect two _SEP separated values per line:
            raise ParseError(
                '{} is not a valid stats file'.format(file_handle.name))
        # extract key and value, then split the key:
        putative_key, putative_val =  fields[0], fields[1]
        key_fields = putative_key.split(_KEY_SEP)
        # we should impose a rigid structure of .stats or redo it:
        if len(key_fields)==1:
            key = key_fields[0]
            if key in stat_from_file:
                stat_from_file[key] = int(fields[1])
            else:
                raise ParseError(
                    '{} is not a valid stats file: unknown field {} detected'.format(file_handle.name,key))
        else:
            # in this case key must be in ['pair_types','chrom_freq','dist_freq','dedup']
            # get the first 'key' and keep the remainders in 'key_fields'
            key = key_fields.pop(0)
            if key in ['pair_types', 'dedup']:
                # assert there is only one element in key_fields left:
                # 'pair_types' and 'dedup' treated the same
                if len(key_fields) == 1:
                    stat_from_file[key][key_fields[0]] = int(fields[1])
                else:
                    raise ParseError(
                        '{} is not a valid stats file: {} section implies 1 identifier'.format(file_handle.name,key))

            elif key == 'chrom_freq':
                # assert remaining key_fields == [chr1, chr2]:
                if len(key_fields) == 2:
                    stat_from_file[key][tuple(key_fields)] = int(fields[1])
                else:
                    raise ParseError(
                        '{} is not a valid stats file: {} section implies 2 identifiers'.format(file_handle.name,key))

            elif key == 'dist_freq':
                # assert that last element of key_fields is the 'directions'
                if len(key_fields) == 2:
                    # assert 'dirs' in ['++','--','+-','-+']
                    dirs = key_fields.pop()
                    # there is only genomic distance range of the bin that's left:
                    bin_range, = key_fields
                    # extract left border of the bin "1000000+" or "1500-6000":
                    dist_bin_left = (bin_range.strip('+')
                        if bin_range.endswith('+')
                        else bin_range.split('-')[0])
                    # get the index of that bin:
                    bin_idx = np.searchsorted(_dist_bins, int(dist_bin_left), 'right') - 1
                    # store corresponding value:
                    stat_from_file[key][dirs][bin_idx] = int(fields[1])
                else:
                    raise ParseError(
                        '{} is not a valid stats file: {} section implies 2 identifiers'.format(file_handle.name,key))
            else:
                raise ParseError(
                    '{} is not a valid stats file: unknown field {} detected'.format(file_handle.name,key))
    # return PairCounter from a non-empty dict:
    #
    # add some fractions:
    stat_from_file['frac_unmapped'] = stat_from_file['total_unmapped']/stat_from_file['total']*100.
    stat_from_file['frac_single_sided_mapped'] = stat_from_file['total_single_sided_mapped']/stat_from_file['total']*100.
    # total_mapped = total_dups + total_nodups
    # should the following be divided by mapped or by total ?!
    stat_from_file['frac_mapped'] = stat_from_file['total_mapped']/stat_from_file['total']*100.
    stat_from_file['frac_dups'] = stat_from_file['total_dups']/stat_from_file['total']*100.
    ########################################
    # the rest of stats are based on nodups:
    ########################################
    stat_from_file['cis_percent'] = stat_from_file['cis']/stat_from_file['total_nodups']*100.

    #
    return stat_from_file



class MultiqcModule(BaseMultiqcModule):
    """This MultiQC module parses various
    stats produced by pairtools."""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='pairtools', anchor='pairtools',
        href="https://github.com/mirnylab/pairtools",
        info="pairtools is a command-line framework for processing sequencing data"
            " generated with Chromatin Conformation Capture based experiments:"
            " pairtools can handle pairs of short-reads aligned to a reference genome,"
            " extract 3C-specific information and perform common tasks, such as sorting,"
            " filtering and deduplication.")

        # Find and load any pairtools stats summaries
        self.pairtools_stats = dict()
        for f in self.find_log_files('pairtools', filehandles=True):
            s_name = f['s_name']
            self.pairtools_stats[s_name] = self.parse_pairtools_stats(f)

        # Filter to strip out ignored sample names
        self.pairtools_stats = self.ignore_samples(self.pairtools_stats)

        if len(self.pairtools_stats) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.pairtools_stats)))
        # self.write_data_file(self.pairtools_stats, 'multiqc_pairtools')

        # max total reads for general stats:
        self.max_total_reads = 0
        for s_name in self.pairtools_stats:
            self.max_total_reads = \
                max(self.max_total_reads, self.pairtools_stats[s_name]['total'])

        self.pairtools_general_stats()

        # Report sections
        self.add_section (
            name = 'Pairs alignment status',
            anchor = 'pair-types',
            description="Number of pairs classified according to their alignment status,"
                        " including uniquely mapped (UU), unmapped (NN), duplicated (DD), and others.",
            helptext = '''For further details check
                        <a href=\"https://pairtools.readthedocs.io/en/latest/formats.html#pair-types\" > pairtools</a>
                        documentation.''',
            plot = self.pair_types_chart()
        )

        self.add_section (
            name = 'Pre-filtered pairs grouped by genomic separations',
            anchor = 'cis-ranges-trans',
            description="Distribution of pre-filtered pairs (UU, UR and RU) by genomic"
                        " separations for <it>cis-</it>pairs and <it>trans-</it>pairs as a separate group.",
            helptext = '''Pre-filtered read pairs might still include artifacts:
            Short-range cis-pairs are typically enriched in technical artifacts, such as self-circles, dangling-ends, etc.
            High fraction of trans interactions typically suggests increased noise levels''',
            plot = self.pairs_by_cisrange_trans()
        )


        self.add_section (
            name = 'Frequency of interactions as a function of genomic separation',
            anchor = 'scalings-plots',
            description="Frequency of interactions (pre-filtered pairs) as a function"
                        " of genomic separation, known as \"scaling plots\", P(s)."
                        " Click on an individual curve to reveal P(s) for different"
                        " read pair orientations.",
            helptext = '''Short-range cis-pairs are typically enriched in technical artifacts.
            Frequency of interactions for read pairs of different orientations
            ++,+-,-+ and -- (FF, FR, RF, RR) provide insight into these technical artifacts.
            Different technical artifacts manifest themselves with only single type of read orientation
            (dangling-ends - FR, self-circles - RF). Thus enrichment of FR/RF pairs at a given genomic
            separation can hint at the level of contamination.''',
            plot = self.pairs_with_genomic_separation()
        )


        self.add_section (
            name = 'Fraction of read pairs by strand orientation',
            anchor = 'read-orientation',
            description="Number of interactions (pre-filtered pairs) reported for every type"
                        " of read pair orientation. Numbers are reported for different"
                        " ranges of genomic separation and combined.",
            helptext = '''Short-range cis-pairs are typically enriched in technical artifacts.
            Frequency of interactions for read pairs of different orientations
            ++,+-,-+ and -- (FF, FR, RF, RR) provide insight into these technical artifacts.
            Different technical artifacts manifest themselves with only single type of read orientation
            (dangling-ends - FR, self-circles - RF). Thus enrichment of FR/RF pairs at a given genomic
            separation can hint at the level of contamination.''',
            plot = self.pairs_by_strand_orientation()
        )


        self.add_section (
            name = 'Pre-filtered pairs grouped by chromosomes',
            anchor = 'pairs-by-chroms',
            description="Number of pre-filtered interactions (pairs) within a single chromosome"
                        " or for a pair of chromosomes.",
            helptext = '''Numbers of pairs are normalized by the total number of pre-filtered pairs per sample.
            Number are reported only for chromosomes/pairs that have >1% of pre-filtered pairs.
            [THERE SEEM TO BE A BUG IN MULTIQC HEATMAP - OBVIOUS WHEN USE HIGHLIGHTING,RENAMING ETC]''',
            plot = self.pairs_by_chrom_pairs()
        )


    def parse_pairtools_stats(self, f):
        """ Parse a pairtools summary stats """
        # s_name = f['s_name']
        # f_name = f['fn']
        # log.info("parsing {} {} ...".format(s_name,f_name))
        f_handle = f['f']
        return _from_file(f_handle)


    def pair_types_chart(self):
        """ Generate the pair types report """

        _report_field = "pair_types"

        known_keys = ['UU', 'RU', 'UR', 'WW', 'DD', 'MR', 'MU', 'MM', 'NM', 'NU', 'NN', 'XX']
        matching_colors = ['#33a02c',
                            '#b2df8a',
                            '#a6cee3',
                            '#1f78b4',
                            '#fff900',
                            '#fb9a99',
                            '#fdbf6f',
                            '#ff7f00',
                            '#e31a1c',
                            '#cab2d6',
                            '#6a3d9a',
                            '#000000']

        # Construct a data structure for the plot
        # and keep track of the observed keys - for nice visuals:
        _data = dict()
        unknown_keys = set()
        for s_name in self.pairtools_stats:
            _data[s_name] = dict()
            for k in self.pairtools_stats[s_name][_report_field]:
                _data[s_name][k] = self.pairtools_stats[s_name][_report_field][k]
                # updated the list of unseen keys
                if k not in known_keys:
                    unknown_keys.add(k)

        # Specify the order of the different possible categories
        keys_anotated = OrderedDict()
        for key, col in zip(known_keys, matching_colors):
            keys_anotated[key] = {'color': col, 'name': key}
        # we have to take care of "unseen" keys like that
        # otherwise they are going to be present on the bar-chart
        for key in unknown_keys:
            keys_anotated[key] = {'name': key}


        # Config for the plot
        config = {
            'id': 'pair_types',
            'title': 'pairtools: pair types report',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        # self.pair_types_chart_data = _data
        return bargraph.plot(_data, keys_anotated, pconfig=config)


    def pairs_by_cisrange_trans(self):
        """ cis-pairs split into several ranges
        of genomic separation, and trans-category."""

        cis_dist_pattern = r"cis_(\d+)kb\+"

        # Construct a data structure for the plot
        _data = dict()
        _tmp = dict()
        for s_name in self.pairtools_stats:
            sample_stats = self.pairtools_stats[s_name]
            _tmp[s_name] = dict()
            _data[s_name] = dict()
            for k in sample_stats:
                # check if a key matches "cis_dist_pattern":
                cis_dist_match = re.fullmatch(cis_dist_pattern, k)
                if cis_dist_match is not None:
                    # extract distance and number of cis-pairs:
                    dist = int(cis_dist_match.group(1))
                    _count = int(sample_stats[k])
                    # fill in _tmp
                    _tmp[s_name][dist] = _count
            # after all cis-dist ones are extracted,
            # undo the cumulative distribution:
            prev_counter = copy(sample_stats['cis'])
            prev_dist = 0
            sorted_keys = []
            for dist in sorted(_tmp[s_name].keys()):
                cis_pair_counter = prev_counter - _tmp[s_name][dist]
                dist_range_key = "cis: {}-{}kb".format(prev_dist,dist)
                _data[s_name][dist_range_key] = cis_pair_counter
                prev_dist = dist
                prev_counter = _tmp[s_name][dist]
                sorted_keys.append(dist_range_key)
            # and the final range max-dist+ :
            cis_pair_counter = _tmp[s_name][dist]
            dist_range_key = "cis: {}+ kb".format(prev_dist)
            _data[s_name][dist_range_key] = prev_counter
            sorted_keys.append(dist_range_key)

            # add trans here as well ...
            dist_range_key = "trans"
            _data[s_name][dist_range_key] = copy(sample_stats['trans'])
            sorted_keys.append(dist_range_key)

        # add color to sorted keys:
        colors_genomic_separations = ['#8c2d04', # short-range cis
                                      '#cc4c02',
                                      '#ec7014',
                                      '#fe9929',
                                      '#fec44f',
                                      '#fee391',
                                      '#ffffd4'] #long-range cis
        key_dict = OrderedDict()
        for col,key in zip_longest(colors_genomic_separations, sorted_keys):
            if col is not None:
                key_dict[key] = {'color': col, 'name': key}
            else:
                # use auto-color when we run out of prescribed ones:
                key_dict[key] = {'name': key}

        # Config for the plot
        config = {
            'id': 'pair_cis_ranges',
            'title': 'pairtools: cis pairs broken into ranges',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        return bargraph.plot(_data, key_dict, pconfig=config)


    def pairs_by_strand_orientation(self):
        """ number of cis-pairs with genomic separation """

        _report_field = "dist_freq"

        dist_edges = [100,500,1000,2000,5000,10000,15000,20000]

        orientation_keys = ['++','-+','+-','--']
        _sorted_keys = ['FF','RF','FR','RR']
        _matching_colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3']

        distances_start = [0]+dist_edges
        distances_end = dist_edges

        dist_ranges = list(zip_longest(distances_start,distances_end))

        # Construct a data structure for the plot
        _data = dict()
        for key in dist_ranges:
            _data[key] = dict()
            for s_name in self.pairtools_stats:
                _data[key][s_name] = dict()


        for s_name in self.pairtools_stats:
            # extract scalings data structure per sample:
            sample_bins = self.pairtools_stats[s_name]['dist_bins']
            sample_dist_freq = self.pairtools_stats[s_name][_report_field]
            # given a set of fixed distance ranges extract FF,RR,FR,RF:
            for start, end in dist_ranges:
                # determine start index, according to "dist_bins"
                start_idx = np.searchsorted(sample_bins, start) if start else 0
                # determine end index, according to "dist_bins"
                end_idx = np.searchsorted(sample_bins, end) if end is not None else None
                # calculate ratios of FF FR RF RR ...
                for orient in orientation_keys:
                    # slice of the scaling for a given range of distances:
                    sliced_data = sample_dist_freq[orient][start_idx:end_idx]
                    _data[(start,end)][s_name][orient] = np.sum(sliced_data.astype(float))

        kb=1000
        data_labels = []
        for start,end in dist_ranges:
            if start == 0:
                data_labels.append(f"<{end}")
            elif end is None:
                # switch to kb if needed
                start_str = f"{start//kb}kb" if start//kb else f"{start}"
                data_labels.append(f">{start_str}")
            else:
                # switch to kb if needed
                if start//kb and end//kb:
                    data_labels.append(f"{start//kb}-{end//kb} kb")
                else:
                    end_str = f"{end//kb}kb" if end//kb else f"{end}"
                    data_labels.append(f"{start}-{end_str}")


        # Config for the plot
        config = {
            'id': 'pair_by_orient_cis_ranges',
            'title': 'pairtools: cis pairs broken into ranges and read orintations',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads',
            'data_labels': data_labels
        }

        # annotate read orientations with nice colors:
        keys_annotated = OrderedDict()
        for key, name, col in zip(orientation_keys, _sorted_keys, _matching_colors):
            keys_annotated[key] = {'color': col, 'name': name}

        return bargraph.plot(
                [_data[_d] for _d in dist_ranges],
                [keys_annotated for _d in dist_ranges],
                pconfig=config)



    # dist_freq/56234133-100000000/-+
    def pairs_with_genomic_separation(self):
        """ number of cis-pairs with genomic separation """

        def _contact_areas(distbins, scaffold_length=2_000_000_000):
            distbins = distbins.astype(float)
            scaffold_length = float(scaffold_length)
            outer_areas = np.maximum(scaffold_length - distbins[:-1], 0) ** 2
            inner_areas = np.maximum(scaffold_length - distbins[1:], 0) ** 2
            return 0.5 * (outer_areas - inner_areas)

        _report_field = "dist_freq"

        # # Construct a data structure for the plot
        # _data_std = dict()
        # _data_spread = dict()
        _data_mean = dict()
        _data_FF = dict()
        _data_FR = dict()
        _data_RR = dict()
        _data_RF = dict()
        for s_name in self.pairtools_stats:
            # _data_std[s_name] = dict()
            # _data_spread[s_name] = dict()
            _data_mean[s_name] = dict()
            _data_FF[s_name] = dict()
            _data_FR[s_name] = dict()
            _data_RR[s_name] = dict()
            _data_RF[s_name] = dict()
            # pre-calculate geom-mean of dist-bins for P(s):
            _dist_bins = self.pairtools_stats[s_name]['dist_bins']
            _areas = _contact_areas( _dist_bins, scaffold_length=2_000_000_000_000 )
            _dist_bins_geom = []
            for i in range(len(_dist_bins)-1):
                geom_dist = np.sqrt(np.prod(_dist_bins[i:i+2]))
                _dist_bins_geom.append(int(geom_dist))
            # _dist_bins_geom are calcualted
            sample_dist_freq = self.pairtools_stats[s_name][_report_field]

            dir_mean = np.sum([
                        sample_dist_freq["++"],
                        sample_dist_freq["+-"],
                        sample_dist_freq["-+"],
                        sample_dist_freq["--"]
                                ],axis=0)[1:].astype(float)
            # ++,--,+-,-+ spread was used before ...

            dir_FF = np.asarray(sample_dist_freq["++"])[1:].astype(float)
            dir_FR = np.asarray(sample_dist_freq["+-"])[1:].astype(float)
            dir_RF = np.asarray(sample_dist_freq["-+"])[1:].astype(float)
            dir_RR = np.asarray(sample_dist_freq["--"])[1:].astype(float)

            dir_mean /= _areas

            dir_FF /= _areas
            dir_FR /= _areas
            dir_RF /= _areas
            dir_RR /= _areas
            #
            # fill in the data ...
            for i, (k,v1) in enumerate(zip(_dist_bins_geom, dir_mean)):
                if i>3:
                    # _data_std[s_name][k] = v1
                    # _data_spread[s_name][k] = v2
                    _data_mean[s_name][k] = v1


            for i, (k,_FF,_FR,_RF,_RR) in enumerate(zip(_dist_bins_geom, dir_FF, dir_FR, dir_RF, dir_RR)):
                if i>3:
                    _data_FF[s_name][k] = _FF
                    _data_FR[s_name][k] = _FR
                    _data_RF[s_name][k] = _RF
                    _data_RR[s_name][k] = _RR

        # # Specify the order of the different possible categories
        # keys = sorted_keys
        # keys['Not_Truncated_Reads'] = { 'color': '#2f7ed8', 'name': 'Not Truncated' }
        # keys['Truncated_Read']      = { 'color': '#0d233a', 'name': 'Truncated' }

        pconfig = {
            'id': 'broom_plot',
            'title': 'Pairs by distance and by read orientation',
            # 'ylab': 'Counts',
            'xlab': 'Genomic separation (bp)',
            'xLog': True,
            'yLog': True,
            # 'xDecimals': False,
            # 'ymin': 0,
            # 'tt_label': '<b>{point.x} bp trimmed</b>: {point.y:.0f}',
            'data_labels': [{'name': 'P(s)', 'ylab': 'frequency of interactions'},
                            {'name': 'FF', 'ylab': 'frequency of interactions'},
                            {'name': 'FR', 'ylab': 'frequency of interactions'},
                            {'name': 'RF', 'ylab': 'frequency of interactions'},
                            {'name': 'RR', 'ylab': 'frequency of interactions'}],
            # 'data_labels': [{'name': 'Spread: std', 'ylab': 'frequency of interactions'},
            #                 {'name': 'Spread: max-min', 'ylab': 'frequency of interactions'},
            #                 {'name': 'P(s): sum', 'ylab': 'frequency of interactions'}],
            'click_func': self.plot_single()
        }

        # plot = linegraph.plot(self.cutadapt_length_counts, pconfig)

        return linegraph.plot([_data_mean, _data_FF, _data_FR, _data_RF, _data_RR], pconfig=pconfig)
        # return linegraph.plot(_data_mean, pconfig=pconfig)


    def plot_single(self):
        """ Return JS code required for plotting a single sample
        RSeQC plot. Attempt to make it look as much like the original as possible.
        Note: this code is injected by `eval(str)`, not <script type="text/javascript"> """

        return """
    function(e){
        // In case of repeated modules: #rseqc_junction_saturation_plot, #rseqc_junction_saturation_plot-1, ..
        var rseqc_junction_saturation_plot = $(e.currentTarget).closest('.hc-plot');
        var rseqc_junction_saturation_plot_id = rseqc_junction_saturation_plot.attr('id');
        var junction_sat_single_hint = rseqc_junction_saturation_plot.closest('.mqc-section').find('#junction_sat_single_hint');


        // Get the three datasets for this sample
        var data = [
            {'name': 'FF'},
            {'name': 'FR'},
            {'name': 'RF'},
            {'name': 'RR'},
        ];
        var k = 0;
        for (var i = 1; i < 5; i++) {
            var ds = mqc_plots[rseqc_junction_saturation_plot_id]['datasets'][i];
            for (k = 0; k < ds.length; k++){
                if(ds[k]['name'] == this.series.name){
                    data[i-1]['data'] = JSON.parse(JSON.stringify(ds[k]['data']));
                    break;
                }
            }
        }

        // Create single plot div, and hide overview
        var newplot = $('<div id="rseqc_junction_saturation_single"> \
            <div id="rseqc_junction_saturation_single_controls"> \
              <button class="btn btn-primary btn-sm" id="rseqc-junction_sat_single_return"> \
                Return to overview \
              </button> \
              <div class="btn-group btn-group-sm"> \
                <button class="btn btn-default rseqc-junction_sat_single_prevnext" data-action="prev">&laquo; Prev</button> \
                <button class="btn btn-default rseqc-junction_sat_single_prevnext" data-action="next">Next &raquo;</button> \
              </div> \
            </div> \
            <div class="hc-plot-wrapper"> \
              <div class="hc-plot hc-line-plot"> \
                <small>loading..</small> \
              </div> \
            </div> \
          </div>');
        var pwrapper = rseqc_junction_saturation_plot.parent().parent();
        newplot.insertAfter(pwrapper).hide().slideDown();
        pwrapper.slideUp();
        junction_sat_single_hint.slideUp();

        // Listener to return to overview
        newplot.find('#rseqc-junction_sat_single_return').click(function(e){
          e.preventDefault();
          newplot.slideUp(function(){
            $(this).remove();
          });
          pwrapper.slideDown();
          junction_sat_single_hint.slideDown();
        });

        // Listeners for previous / next plot
        newplot.find('.rseqc-junction_sat_single_prevnext').click(function(e){
          e.preventDefault();
          if($(this).data('action') == 'prev'){
            k--;
            if(k < 0){
              k = mqc_plots[rseqc_junction_saturation_plot_id]['datasets'][1].length - 1;
            }
          } else {
            k++;
            if(k >= mqc_plots[rseqc_junction_saturation_plot_id]['datasets'][1].length){
              k = 0;
            }
          }
          var hc = newplot.find('.hc-plot').highcharts();
          for (var i = 1; i < 5; i++) {
              hc.series[i-1].setData(mqc_plots[rseqc_junction_saturation_plot_id]['datasets'][i][k]['data'], false);
          }
          var ptitle = 'frequency of interactions: '+mqc_plots[rseqc_junction_saturation_plot_id]['datasets'][1][k]['name'];
          hc.setTitle({text: ptitle});
          hc.redraw({ duration: 200 });
        });

        // Plot the single data
        newplot.find('.hc-plot').highcharts({
          chart: {
            type: 'line',
            zoomType: 'x'
          },
          colors: ['red','green','blue','orange'],
          title: {
            text: 'frequency of interaction for '+this.series.name,
            x: 30 // fudge to center over plot area rather than whole plot
          },
          xAxis: {
            title: { text: 'Genomic separation (bp)' },
            allowDecimals: false,
            type: 'logarithmic',
          },
          yAxis: {
            title: { text: 'frequency of interactions' },
            type: 'logarithmic',
          },
          legend: {
            floating: true,
            layout: 'vertical',
            align: 'left',
            verticalAlign: 'top',
            x: 60,
            y: 40
          },
          tooltip: {
            shared: true,
            crosshairs: true,
            headerFormat: '<strong>genomic separation: {point.key} bp</strong><br/>'
          },
          plotOptions: {
            series: {
              animation: false,
              lineWidth: 2,
              marker: {
                lineColor: null,
                fillColor: 'transparent',
                lineWidth: 0,
                symbol: 'circle'
              },
            }
          },
          exporting: { buttons: { contextButton: {
            menuItems: window.HCDefaults.exporting.buttons.contextButton.menuItems,
            onclick: window.HCDefaults.exporting.buttons.contextButton.onclick
          } } },
          series: data
        });
    }
    """



    # chrom_freq/chr1/chrX ...
    def pairs_by_chrom_pairs(self):
        """ number of pairs by chromosome pairs """

        _report_field = "chrom_freq"

        # figure infer list of chromosomes (beware of scaffolds):
        # tuple(key_fields)
        _chromset = set()
        for s_name in self.pairtools_stats:
            _chrom_freq_sample = \
                self.pairtools_stats[s_name][_report_field]
            # unzip list of tuples:
            _chroms1, _chroms2 = list(
                    zip(*_chrom_freq_sample.keys())
                )
            _chromset |= set(_chroms1)
            _chromset |= set(_chroms2)
        # done:
        _chroms = sorted(list(_chromset))

        # Construct a data structure for the plot
        _data = dict()
        for s_name in self.pairtools_stats:
            _data[s_name] = dict()
            _chrom_freq_sample = \
                self.pairtools_stats[s_name][_report_field]
            # go over chroms:
            for c1,c2 in combinations_with_replacement( _chroms, 2):
                # record the chromosome combination:
                _chrom_combo = (c1,c2)
                # _num_pairs calculations:
                if (c1,c2) in _chrom_freq_sample:
                    _num_pairs = _chrom_freq_sample[(c1,c2)]
                elif (c2,c1) in _chrom_freq_sample:
                    _num_pairs = _chrom_freq_sample[(c2,c1)]
                    # _chrom_combo = (c2,c1)
                else:
                    _num_pairs = 0
                # let's filter by # of pairs - by doing some masking ...
                if _num_pairs < 0.0007*self.pairtools_stats[s_name]['cis']:
                    _num_pairs = 0
                else:
                    # we'll try to normalize it afterwards ...
                    _num_pairs /= (self.pairtools_stats[s_name]['cis']+self.pairtools_stats[s_name]['trans'])
                    # pass
                _data[s_name][_chrom_combo] = _num_pairs

        # now we need to filter 0 cells ...
        # prepare for the heatmap:
        xcats = sorted([ (c1, c2) for c1, c2 in combinations_with_replacement( _chroms, 2) ])
        xcats_names = [f"{c1}-{c2}" for c1,c2 in xcats]
        # try samples as x-category ...
        ycats = sorted(_data)
        the_data = [ [ _data[s][k] for k in xcats ] for s in ycats ]


        # check if there are any zeros in the column (i.e. for a given chrom pair) ...
        mask = np.all(the_data, axis=0)
        the_data_filt = np.asarray(the_data)[:,mask]
        # # mean over columns to sort ...
        sorted_idx = the_data_filt.mean(axis=0).argsort()
        return heatmap.plot(
                the_data_filt[:,sorted_idx].tolist(),
                np.array(xcats_names)[mask][sorted_idx].tolist(),
                ycats)#, pconfig)


    def pairtools_general_stats(self):
        """ Add columns to General Statistics table """
        headers = OrderedDict()
        headers['total'] = {
            'title': '{} read pairs'.format(config.read_count_prefix),
            'description': 'Total read pairs ({})'.format(config.read_count_desc),
            'min': 0,
            'modify': lambda x: x * config.read_count_multiplier,
            'scale': 'Blues',
        }
        headers['frac_unmapped'] = {
            'title': '% unmapped',
            'description': '% of pairs (w.r.t. total) with both sides unmapped',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'OrRd',
        }
        headers['frac_single_sided_mapped'] = {
            'title': '% single-side mapped',
            'description': '% of pairs (w.r.t. total) with one side mapped',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn',
        }
        headers['frac_mapped'] = {
            'title': '% both-side mapped',
            'description': '% of pairs (w.r.t. total) with both sides mapped',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn',
        }
        headers['frac_dups'] = {
            'title': '% duplicated',
            'description': '% of duplicated pairs (w.r.t. mapped)',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'OrRd',
        }
        headers['cis_percent'] = {
            'title': '% cis',
            'description': '% of cis-pairs (w.r.t mapped)',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn',
        }
        self.general_stats_addcols(self.pairtools_stats, headers, 'pairtools')
