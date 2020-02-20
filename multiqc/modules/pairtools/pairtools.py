#!/usr/bin/env python

""" MultiQC module to parse stats output from pairtools """

import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule


from itertools import combinations, combinations_with_replacement
import re
import numpy as np
from copy import copy, deepcopy
from multiqc.plots import bargraph, linegraph, heatmap

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
    stat_from_file['frac_nodups'] = stat_from_file['total_nodups']/stat_from_file['total']*100.
    ########################################
    # the rest of stats are based on nodups:
    ########################################
    stat_from_file['cis_percent'] = stat_from_file['cis']/stat_from_file['total_nodups']*100.
    stat_from_file['trans_percent'] = stat_from_file['trans']/stat_from_file['total_nodups']*100.

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

        self.pairtools_general_stats()

        # Report sections
        self.add_section (
            name = 'Pairs alignment status report',
            description="Number of pairs clasified according to their alignment status,"
                        " including uniquely mapped (UU), unmapped (NN), duplicated (DD), and others."
                        " For further details check"
                        " <a href=\"https://pairtools.readthedocs.io/en/latest/formats.html#pair-types\" > pairtools</a>"
                        " documentation.",
            anchor = 'pair-types',
            plot = self.pair_types_chart()
        )

        self.add_section (
            name = 'Usable pairs grouped by genomic separations',
            anchor = 'cis-ranges',
            description="Distribution of usable pairs (UU, UR and RU) by their genomic"
                        " separations for cis-pairs and trans-pairs as a separate group."
                        " Short-range cis-pairs are typically enriched in technical artifacts.",
            plot = self.cis_breakdown_chart()
        )

        self.add_section (
            name = 'Frequency of interactions as a function of genomic separation',
            anchor = 'pairs dirs',
            description="Frequency of interactions (usable pairs) reported"
                        " as a function of genomic separation, known as \"scaling plots\", P(s)."
                        " Frequency of interactions for pairs of different strand orientations"
                        " ++,+-,-+ and -- (FF, FR, RF, RR) provide insight into technical artifacts:"
                        " typically manifested as unequal fractions of FF, FR, RF, RR at a given separation."
                        " In order to avoid showing 4 lines per sample we plot std(FF,FR,RF,RR), as a measure"
                        " of \"divergence\" between pairs of different strand orientations."
                        " [THIS SECTION IS MY BIGGEST STRUGGLE - NOT SURE IF STANDARD DEVIATION DOES A GOOD JOB"
                        " AT SHOWING THE DIVERGENCE OF A \"BROOM\"-PLOT, CONSIDERED MAX(FF,FR,RF,RR)-MIN(FF,FR,RF,RR)"
                        " BEST THING WOULD BE TO SHOW 4-LINES FF,FR,RF,RR, AFTER CLICKING ON A GIVEN SAMPLE"
                        " - WOULD REQUIRE SOME JAVASCRIPTING - PROBABLY NOT FOR THE FIRST ITERATION... ]",
            plot = self.pairs_with_genomic_separation()
        )

        self.add_section (
            name = 'Usable pairs grouped by chromosomes',
            anchor = 'pairs chrom/chroms ...',
            description="Number of usable pairs interacting within each chromosome"
                        " or for a combination of chromosomes."
                        " Numbers of pairs are normalized by the total number of"
                        " usable pairs per sample."
                        " Number are reported only for chromosomes/combinations that have >1% of usable pairs."
                        " [THERE SEEM TO BE A BUG IN MULTIQC HEATMAP - OBVIOUS WHEN USE HIGHLIGHTING,RENAMING ETC]",
            plot = self.pairs_by_chrom_pairs()
        )

        self.add_section (
            name = 'Inter-sample correlation [ATTEMPT]',
            anchor = 'pairs chrom/chroms corr...',
            description="An attempt to characterize inter-sample similarity"
                        " using different available metrics:"
                        " distribution of pairs per chromosome pair"
                        " Samples with translocations would 'cluster' away from the 'normal' ones with this metric..."
                        " If we add chromsizes as part of stats, we could display # of cis-pairs ~ chromlen, or"
                        " # of trans-pairs ~inter-chromosomal 'area' (these are supposed to be ~linear relationships)"
                        " so, samples with translocations and other stuff would pop-up ....(maybe)"
                        " Another thing to do (bassed on this) - show 'coverage' per chromosome ....",
            plot = self.samples_similarity_chrom_pairs()
        )

        self.add_section (
            name = 'Inter-sample correlation [ATTEMPT#2]',
            anchor = 'pairs by separation groups...',
            description="An attempt #2 to characterize inter-sample similarity"
                        " using different available metrics:"
                        " distribution of pairs by distance (and trans)"
                        " Samples with translocations would 'cluster' away from the 'normal' ones with this metric..."
                        " Overall this should show similar 'clustering' as in attempt #1 ...",
            plot = self.samples_similarity_cis_pairs()
        )


        self.add_section (
            name = 'Inter-sample correlation [ATTEMPT#3]',
            anchor = 'pairs by alignment status...',
            description="An attempt #3 to characterize inter-sample similarity"
                        " using different available metrics:"
                        " distribution of pairs alignment status (UU,NN,MM,...)",
            plot = self.samples_similarity_mapped_pairs()
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

        # Construct a data structure for the plot
        _data = dict()
        for s_name in self.pairtools_stats:
            _data[s_name] = dict()
            for k in self.pairtools_stats[s_name][_report_field]:
                _data[s_name][k] = self.pairtools_stats[s_name][_report_field][k]

        # # Specify the order of the different possible categories
        # keys = OrderedDict()
        # keys['Not_Truncated_Reads'] = { 'color': '#2f7ed8', 'name': 'Not Truncated' }
        # keys['Truncated_Read']      = { 'color': '#0d233a', 'name': 'Truncated' }

        # Config for the plot
        config = {
            'id': 'pair_types',
            'title': 'pairtools: pair types report',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        self.pair_types_chart_data = _data

        return bargraph.plot(_data, pconfig=config)


    def cis_breakdown_chart(self):
        """ cis-pairs split into several ranges """

        cis_dist_pattern = r"cis_(\d+)kb\+"

        # Construct a data structure for the plot
        _data = dict()
        _datawtrans = dict()
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
            # after all cis-dist ones are extracted:
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
            _datawtrans[s_name] = dict(_data[s_name])
            dist_range_key = "trans"
            _datawtrans[s_name][dist_range_key] = copy(sample_stats['trans'])
            sorted_keys.append(dist_range_key)

        # # Specify the order of the different possible categories
        # keys = sorted_keys
        # keys['Not_Truncated_Reads'] = { 'color': '#2f7ed8', 'name': 'Not Truncated' }
        # keys['Truncated_Read']      = { 'color': '#0d233a', 'name': 'Truncated' }

        self.pairs_breakdown = _datawtrans

        # Config for the plot
        config = {
            'id': 'pair_cis_ranges',
            'title': 'pairtools: cis pairs broken into ranges',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads',
            # 'logswitch': True - useless for that
            # 'data_labels': ['cis-only','cis-n-trans']
        }

        return bargraph.plot(_datawtrans, sorted_keys, pconfig=config)


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

        # Construct a data structure for the plot
        _data_std = dict()
        _data_spread = dict()
        _data_mean = dict()
        for s_name in self.pairtools_stats:
            _data_std[s_name] = dict()
            _data_spread[s_name] = dict()
            _data_mean[s_name] = dict()
            # pre-calculate geom-mean of dist-bins for P(s):
            _dist_bins = self.pairtools_stats[s_name]['dist_bins']
            _areas = _contact_areas( _dist_bins, scaffold_length=2_000_000_000_000 )
            _dist_bins_geom = []
            for i in range(len(_dist_bins)-1):
                geom_dist = np.sqrt(np.prod(_dist_bins[i:i+2]))
                _dist_bins_geom.append(int(geom_dist))
            # _dist_bins_geom are calcualted
            sample_dist_freq = self.pairtools_stats[s_name][_report_field]
            # that's how it is for now ...
            # stat_from_file[key][dirs][bin_idx] = int(fields[1])
            dir_std = np.std([
                        sample_dist_freq["++"],
                        sample_dist_freq["+-"],
                        sample_dist_freq["-+"],
                        sample_dist_freq["--"]
                                ],axis=0)[1:].astype(float)
            # consider using max-min instead ...
            # amount of information sum(f*log(f) ) ?!...

            # data spread max(++,+-,-+,--) - min(++,+-,-+,--)
            dir_spread = \
                    np.max([
                        sample_dist_freq["++"],
                        sample_dist_freq["+-"],
                        sample_dist_freq["-+"],
                        sample_dist_freq["--"]
                                ],axis=0)[1:].astype(float) - \
                    np.min([
                        sample_dist_freq["++"],
                        sample_dist_freq["+-"],
                        sample_dist_freq["-+"],
                        sample_dist_freq["--"]
                                ],axis=0)[1:].astype(float)

            dir_mean = np.sum([
                        sample_dist_freq["++"],
                        sample_dist_freq["+-"],
                        sample_dist_freq["-+"],
                        sample_dist_freq["--"]
                                ],axis=0)[1:].astype(float)
            # / self.pairtools_stats[s_name]["cis_1kb+"]

            # dir_std /= _areas#+0.01
            # dir_spread /= _areas#+0.01
            dir_mean /= _areas#+0.01
            #
            # fill in the data ...
            for i,(k,v1,v2,v3) in enumerate(zip(_dist_bins_geom, dir_std, dir_spread, dir_mean)):
                if i>3:
                    _data_std[s_name][k] = v1
                    _data_spread[s_name][k] = v2
                    _data_mean[s_name][k] = v3

        # # Specify the order of the different possible categories
        # keys = sorted_keys
        # keys['Not_Truncated_Reads'] = { 'color': '#2f7ed8', 'name': 'Not Truncated' }
        # keys['Truncated_Read']      = { 'color': '#0d233a', 'name': 'Truncated' }

        pconfig = {
            'id': 'broom plot',
            'title': 'Pairs by distance and by read orientation',
            # 'ylab': 'Counts',
            'xlab': 'Genomic separation (bp)',
            'xLog': True,
            'yLog': True,
            # 'xDecimals': False,
            # 'ymin': 0,
            # 'tt_label': '<b>{point.x} bp trimmed</b>: {point.y:.0f}',
            'data_labels': [{'name': 'Spread: std', 'ylab': 'frequency of interactions'},
                            {'name': 'Spread: max-min', 'ylab': 'frequency of interactions'},
                            {'name': 'P(s): sum', 'ylab': 'frequency of interactions'}]
        }

        # plot = linegraph.plot(self.cutadapt_length_counts, pconfig)

        return linegraph.plot([_data_std, _data_spread, _data_mean], pconfig=pconfig)


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





    # chrom_freq/chr1/chrX ...
    def samples_similarity_chrom_pairs(self):
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

        # cis-only for now:
        # Construct a data structure for the plot
        _data = dict()
        for s_name in self.pairtools_stats:
            _data[s_name] = []
            _chrom_freq_sample = \
                self.pairtools_stats[s_name][_report_field]
            # go over chroms:
            for c1,c2 in combinations_with_replacement( _chroms, 2):
                if (c1,c2) in _chrom_freq_sample:
                    _num_pairs = _chrom_freq_sample[(c1,c2)]
                elif (c2,c1) in _chrom_freq_sample:
                    _num_pairs = _chrom_freq_sample[(c2,c1)]
                else:
                    _num_pairs = 0
                # let's take ALL data in ...
                # we'll try to normalize it afterwards ...
                _num_pairs /= (self.pairtools_stats[s_name]['cis']+self.pairtools_stats[s_name]['trans'])
                # pass
                _data[s_name].append(_num_pairs)

        #


        # prepare for the heatmap:
        ycats = sorted(_data)
        the_data = np.asarray([ _data[_] for _ in ycats ])
        corrs = np.corrcoef(the_data)
        #


        # let's try other metric here as well ...
        the_data2 = [[self.pairs_breakdown[s][k] for k in sorted(self.pairs_breakdown[s])] for s in ycats]
        corrs2 = np.corrcoef(the_data2)


        return heatmap.plot(
                corrs.tolist(),
                ycats,
                ycats)#, pconfig)


    def samples_similarity_cis_pairs(self):
        """ number of pairs by chromosome pairs """

        # prepare for the heatmap:
        ycats = sorted(self.pairs_breakdown)

        # let's try other metric here as well ...
        # normalize the data - i.e. use percentage instead of the counts ...
        the_data = [
                        [
                          float(self.pairs_breakdown[s][k])/sum(self.pairs_breakdown[s].values()) \
                            for k in sorted(self.pairs_breakdown[s])
                        ] \
                      for s in ycats]

        corrs = np.corrcoef(the_data)


        return heatmap.plot(
                corrs.tolist(),
                ycats,
                ycats)#, pconfig)


    def samples_similarity_mapped_pairs(self):
        """ number of pairs by chromosome pairs """

        # prepare for the heatmap:
        ycats = sorted(self.pair_types_chart_data)

        # let's try other metric here as well ...
        # normalize the data - i.e. use percentage instead of the counts ...
        the_data = [
                    [
                     float(self.pair_types_chart_data[s][k])/sum(self.pair_types_chart_data[s].values()) \
                       for k in sorted(self.pair_types_chart_data[s])
                    ] for s in ycats
                   ]
        corrs = np.corrcoef(the_data)


        return heatmap.plot(
                corrs.tolist(),
                ycats,
                ycats)#, pconfig)


    def pairtools_general_stats(self):
        """ Add columns to General Statistics table """
        headers = OrderedDict()
        headers['total'] = {
            'title': 'total',
            'description': 'total number of pairs per sample',
            'min': 0,
        }
        headers['frac_unmapped'] = {
            'title': 'unmapped',
            'description': 'fraction of unmapped',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn-rev',
        }
        headers['frac_single_sided_mapped'] = {
            'title': 'single-side mapped',
            'description': 'fraction of single-side mapped',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn-rev',
        }
        headers['frac_mapped'] = {
            'title': 'both-side mapped',
            'description': 'fraction of both-side mapped',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn-rev',
        }
        headers['frac_dups'] = {
            'title': 'duplicates',
            'description': 'fraction of duplicates',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn-rev',
        }
        headers['frac_nodups'] = {
            'title': 'unique',
            'description': 'fraction of unique',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn-rev',
        }
        headers['cis_percent'] = {
            'title': 'cis',
            'description': 'cis percent',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn-rev',
        }
        # this is redundant
        # headers['trans_percent'] = {
        #     'title': 'trans',
        #     'description': 'trans percent',
        #     'max': 100,
        #     'min': 0,
        #     'suffix': '%',
        #     'scale': 'YlGn-rev',
        # }
        self.general_stats_addcols(self.pairtools_stats, headers, 'pairtools')
        # self.general_stats_addcols(self.pairtools_stats)

