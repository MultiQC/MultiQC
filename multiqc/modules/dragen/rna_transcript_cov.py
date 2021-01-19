#!/usr/bin/env python
from __future__ import print_function

# Initialise the logger
import logging
import re
from collections import defaultdict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


class DragenRnaTranscriptCoverage(BaseMultiqcModule):
    def add_rna_transcript_coverage(self):
        data_by_sample = defaultdict(dict)

        for f in self.find_log_files('dragen/rna_transcript_cov'):
            data = parse_rna_transcript_cov(f)
            if f['s_name'] in data_by_sample:
                log.debug('Duplicate sample name found! Overwriting: {}'.format(f['s_name']))
            self.add_data_source(f, section='stats')
            data_by_sample[f['s_name']] = data

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)

        if not data_by_sample:
            return set()

        self.add_section(
            name='RNA Transcript Coverage',
            anchor='rna-transcript-cov',
            description="""
            RNA transcript coverage.  This is the average coverage at the position along the transcripts.
            """,
            plot=linegraph.plot(data_by_sample, pconfig={
                'id': 'dragen_rna_transcript_cov',
                'title': 'Dragen: RNA Transcript Coverage',
                'ylab': 'Average coverage',
                'xlab': 'Transcript position',
                'categories': True,
                'tt_label': '<b>{point.x}</b>: {point.y:.1f}x',
            })
        )

        return data_by_sample.keys()


def parse_rna_transcript_cov(f):
    """
    *.quant.transcript_coverage.txt
    Percentile	Coverage
    1	0.0518503
    2	0.125753
    3	0.194324
    4	0.254819
    5	0.308795
    """

    f['s_name'] = re.search(r'(.*).quant.transcript_coverage.txt', f['fn']).group(1)

    data = {}
    for line in f['f'].splitlines()[1:]:
        percentile, coverage = line.split()

        try:
            percentile = float(percentile)
        except ValueError:
            pass
        try:
            coverage = float(coverage)
        except ValueError:
            pass
        data[percentile] = coverage

    return data
