#!/usr/bin/env python

""" MultiQC module to parse stats output from pairtools """

import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule

import numpy as np

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
        info="pairtools is a command-line framework to process sequencing data from a Hi-C experiment.")


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

        # # Report sections
        # self.add_section (
        #     name = 'Read Truncation',
        #     anchor = 'hicup-truncating',
        #     plot = self.hicup_truncating_chart()
        # )



    def parse_pairtools_stats(self, f):
        """ Parse a pairtools summary stats """
        #
        # categories that we need to read from stats file
        # ['chrom_freq',
        #  'cis',
        #  'cis_10kb+',
        #  'cis_1kb+',
        #  'cis_20kb+',
        #  'cis_2kb+',
        #  'cis_40kb+',
        #  'cis_4kb+',
        #  'dist_freq',
        #  'pair_types',
        #  'total',
        #  'total_dups',
        #  'total_mapped',
        #  'total_nodups',
        #  'total_single_sided_mapped',
        #  'total_unmapped',
        #  'trans']
        #
        #
        # s_name = f['s_name']
        # f_name = f['fn']
        # log.info("parsing {} {} ...".format(s_name,f_name))
        #
        # just testing displying random fields from stats ...
        f_handle = f['f']
        # _data = dict()
        return _from_file(f_handle)
        # _i = 0
        # num_fields = 8
        # for line in f_handle:
        #     _1,_2 = line.rstrip().split('\t')
        #     _data[_1] = int(_2)
        #     _i += 1
        #     if _i >=num_fields:
        #         break
        # return _data



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


    # def hicup_truncating_chart (self):
    #     """ Generate the HiCUP Truncated reads plot """

    #     # Specify the order of the different possible categories
    #     keys = OrderedDict()
    #     keys['Not_Truncated_Reads'] = { 'color': '#2f7ed8', 'name': 'Not Truncated' }
    #     keys['Truncated_Read']      = { 'color': '#0d233a', 'name': 'Truncated' }

    #     # Construct a data structure for the plot - duplicate the samples for read 1 and read 2
    #     data = {}
    #     for s_name in self.hicup_data:
    #         data['{} Read 1'.format(s_name)] = {}
    #         data['{} Read 2'.format(s_name)] = {}
    #         data['{} Read 1'.format(s_name)]['Not_Truncated_Reads'] = self.hicup_data[s_name]['Not_Truncated_Reads_1']
    #         data['{} Read 2'.format(s_name)]['Not_Truncated_Reads'] = self.hicup_data[s_name]['Not_Truncated_Reads_2']
    #         data['{} Read 1'.format(s_name)]['Truncated_Read'] = self.hicup_data[s_name]['Truncated_Read_1']
    #         data['{} Read 2'.format(s_name)]['Truncated_Read'] = self.hicup_data[s_name]['Truncated_Read_2']

    #     # Config for the plot
    #     config = {
    #         'id': 'hicup_truncated_reads_plot',
    #         'title': 'HiCUP: Truncated Reads',
    #         'ylab': '# Reads',
    #         'cpswitch_counts_label': 'Number of Reads'
    #     }

    #     return bargraph.plot(data, keys, config)

