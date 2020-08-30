
import numpy as np
from collections import OrderedDict


_SEP = '\t'
_KEY_SEP = '/'


# some good-looking presets for pair-type report section:
pairtypes_common = OrderedDict([('UU','#33a02c'),
    ('RU','#b2df8a'),
    ('UR','#a6cee3'),
    ('WW','#1f78b4'),
    ('DD','#fff900'),
    ('MR','#fb9a99'),
    ('MU','#fdbf6f'),
    ('MM','#ff7f00'),
    ('NM','#e31a1c'),
    ('NU','#cab2d6'),
    ('NN','#6a3d9a'),
    ('XX','#000000')])


# add color to sorted keys:
cis_range_colors = ['#8c2d04', # short-range cis
    '#cc4c02',
    '#ec7014',
    '#fe9929',
    '#fec44f',
    '#fee391',
    '#ffffd4'] #long-range cis


# this should be based on chromosome sizes, not an arbitrary number
def contact_areas(distbins, scaffold_length=2_000_000_000):
    distbins = distbins.astype(float)
    scaffold_length = float(scaffold_length)
    outer_areas = np.maximum(scaffold_length - distbins[:-1], 0) ** 2
    inner_areas = np.maximum(scaffold_length - distbins[1:], 0) ** 2
    return 0.5 * (outer_areas - inner_areas)



class ParseError(Exception):
    pass


def read_pairs_stats(file_handle):
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

    # we probably do not need to define these here,
    # just read the 1-key entries and then parse those out:
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
