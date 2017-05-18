#!/usr/bin/env python

""" MultiQC submodule to parse output from Bamtools bam_stat.py
http://bamtools.sourceforge.net/#bam-stat-py """

from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import beeswarm

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find bamtools stats reports and parse their data """

    # Set up vars
    self.bamtools_stats_data = dict()
    regexes = {
        'total_reads': r"Total reads:\s*(\d+)",
        'mapped_reads': r"Mapped reads:\s*(\d+)",
        'mapped_reads_pct': r"Mapped reads:\s*\d+\s+\(([\d\.]+)%\)",
        'forward_strand': r"Forward strand:\s*(\d+)",
        'forward_strand_pct': r"Forward strand:\s*\d+\s+\(([\d\.]+)%\)",
        'reverse_strand': r"Reverse strand:\s*(\d+)",
        'reverse_strand_pct': r"Reverse strand:\s*\d+\s+\(([\d\.]+)%\)",
        'failed_qc': r"Failed QC:\s*(\d+)",
        'failed_qc_pct': r"Failed QC:\s*\d+\s+\(([\d\.]+)%\)",
        'duplicates': r"Duplicates:\s*(\d+)",
        'duplicates_pct': r"Duplicates:\s*\d+\s+\(([\d\.]+)%\)",
        'paired_end': r"Paired-end reads:\s*(\d+)",
        'paired_end_pct': r"Paired-end reads:\s*\d+\s+\(([\d\.]+)%\)",
        'proper_pairs': r"'Proper-pairs'\s*(\d+)",
        'proper_pairs_pct': r"'Proper-pairs'\s*\d+\s+\(([\d\.]+)%\)",
        'both_mapped': r"Both pairs mapped:\s*(\d+)",
        'both_mapped_pct': r"Both pairs mapped:\s*\d+\s+\(([\d\.]+)%\)",
        'read_1': r"Read 1:\s*(\d+)",
        'read_2': r"Read 2:\s*(\d+)",
        'singletons': r"Singletons:\s*(\d+)",
        'singletons_pct': r"Singletons:\s*\d+\s+\(([\d\.]+)%\)",
    }

    # Go through files and parse data using regexes
    for f in self.find_log_files('bamtools/stats'):
        d = dict()
        for k, r in regexes.items():
            r_search = re.search(r, f['f'], re.MULTILINE)
            if r_search:
                d[k] = float(r_search.group(1))

        if len(d) > 0:
            if f['s_name'] in self.bamtools_stats_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            self.add_data_source(f, section='stats')
            self.bamtools_stats_data[f['s_name']] = d

    # Filter to strip out ignored sample names
    self.bamtools_stats_data = self.ignore_samples(self.bamtools_stats_data)

    if len(self.bamtools_stats_data) > 0:

        # Write to file
        self.write_data_file(self.bamtools_stats_data, 'multiqc_bamtools_stats')

        # Add to general stats table
        self.general_stats_headers['duplicates_pct'] = {
            'title': '% Duplicates',
            'description': '% Duplicate Reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'OrRd'
        }
        self.general_stats_headers['mapped_reads_pct'] = {
            'title': '% Mapped',
            'description': '% Mapped Reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn'
        }
        for s_name in self.bamtools_stats_data:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update( self.bamtools_stats_data[s_name] )

        # Make dot plot of counts
        keys = OrderedDict()
        defaults = {
            'min': 0,
            'max': 100,
            'decimalPlaces': 2,
            'suffix': '%'
        }
        num_defaults = {
            'min': 0,
            'modify': lambda x: float(x) / 1000000.0,
            'decimalPlaces': 2
        }

        keys['total_reads'] = dict(num_defaults, **{'title': 'Total reads', 'description': 'Total reads (millions)' });
        keys['mapped_reads_pct'] = dict(defaults, **{'title': 'Mapped reads' })
        keys['forward_strand_pct'] = dict(defaults, **{'title': 'Forward strand' })
        keys['reverse_strand_pct'] = dict(defaults, **{'title': 'Reverse strand' })
        keys['failed_qc_pct'] = dict(defaults, **{'title': 'Failed QC' })
        keys['duplicates_pct'] = dict(defaults, **{'title': 'Duplicates' })
        keys['paired_end_pct'] = dict(defaults, **{'title': 'Paired-end', 'description': 'Paired-end reads' })
        keys['proper_pairs_pct'] = dict(defaults, **{'title': 'Proper-pairs' })
        keys['both_mapped_pct'] = dict(defaults, **{'title': 'Both mapped', 'description': 'Both pairs mapped' })
        keys['read_1'] = dict(num_defaults, **{'title': 'Read 1', 'description': 'Read 1 (millions)' });
        keys['read_2'] = dict(num_defaults, **{'title': 'Read 2', 'description': 'Read 2 (millions)' });
        keys['singletons_pct'] = dict(defaults, **{'title': 'Singletons' })

        self.add_section (
            name = 'Bamtools Stats',
            anchor = 'bamtools-stats',
            plot = beeswarm.plot(self.bamtools_stats_data, keys)
        )

    # Return number of samples found
    return len(self.bamtools_stats_data)

