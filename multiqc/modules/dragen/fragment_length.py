#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict, defaultdict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

# Initialise the logger
import logging

log = logging.getLogger(__name__)


class DragenFragmentLength(BaseMultiqcModule):
    def parse_fragment_length_hist(self):
        """
        #Sample: N_SRR7890889
        FragmentLength,Count
        36,1
        37,0
        38,0
        39,0
        40,0
        41,1
        ...
        39203,0
        39204,0
        39205,1
        #Sample: T_SRR7890936_50pc
        FragmentLength,Count
        53,2
        54,0
        ...
        39316,0
        39317,1
        """

        all_data_by_sample = dict()

        for f in self.find_log_files('dragen/fragment_length_hist'):
            data_by_sample = parse_fragment_length_hist_file(f)
            if data_by_sample:
                for sn, data in data_by_sample.items():
                    if sn in all_data_by_sample:
                        log.debug('Duplicate sample name found! Overwriting: {}'.format(f['s_name']))
                    self.add_data_source(f, section='stats')

                    all_data_by_sample[sn] = data

        # Filter to strip out ignored sample names:
        all_data_by_sample = self.ignore_samples(all_data_by_sample)

        if not all_data_by_sample:
            return
        log.info('Found fragment length histogram for {} samples'.format(len(all_data_by_sample)))

        self.add_section(
            name='Fragment length histogram',
            anchor='dragen-fragment-length-histogram',
            description='Distribution of estimated fragment lengths of mapped reads.',
            plot=linegraph.plot(all_data_by_sample, {
                'id': 'dragen_fragment_length',
                'title': 'Dragen: fragment length histogram',
                'ylab': 'Fraction of reads',
                'xlab': 'Fragment length (bp)',
                'ymin': 0,
                'xmin': 0,
                'tt_label': '<b>{point.x} bp</b>: {point.y}',
            })
        )


def parse_fragment_length_hist_file(f):
    data_by_sample = defaultdict(dict)

    sample = None
    for line in f['f'].splitlines():
        if line.startswith('#Sample'):
            sample = line.split('#Sample: ')[1]
        else:
            assert sample is not None
            frag_len, cnt = line.split(',')
            try:
                frag_len = int(frag_len)
                cnt = int(cnt)
            except ValueError:
                assert line == 'FragmentLength,Count', line
            else:
                data_by_sample[sample][frag_len] = cnt

    return data_by_sample


