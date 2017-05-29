#!/usr/bin/env python

""" MultiQC submodule to parse output from RSeQC bam_stat.py
http://rseqc.sourceforge.net/#bam-stat-py """

from collections import OrderedDict
import logging
import re

from multiqc.plots import beeswarm

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find RSeQC bam_stat reports and parse their data """

    # Set up vars
    self.bam_stat_data = dict()
    regexes = {
        'total_records': r"Total records:\s*(\d+)",
        'qc_failed': r"QC failed:\s*(\d+)",
        'optical_pcr_duplicate': r"Optical/PCR duplicate:\s*(\d+)",
        'non_primary_hits': r"Non primary hits\s*(\d+)",
        'unmapped_reads': r"Unmapped reads:\s*(\d+)",
        'mapq_lt_mapq_cut_non-unique': r"mapq < mapq_cut \(non-unique\):\s*(\d+)",
        'mapq_gte_mapq_cut_unique': r"mapq >= mapq_cut \(unique\):\s*(\d+)",
        'read_1': r"Read-1:\s*(\d+)",
        'read_2': r"Read-2:\s*(\d+)",
        'reads_map_to_sense': r"Reads map to '\+':\s*(\d+)",
        'reads_map_to_antisense': r"Reads map to '-':\s*(\d+)",
        'non-splice_reads': r"Non-splice reads:\s*(\d+)",
        'splice_reads': r"Splice reads:\s*(\d+)",
        'reads_mapped_in_proper_pairs': r"Reads mapped in proper pairs:\s*(\d+)",
        'proper-paired_reads_map_to_different_chrom': r"Proper-paired reads map to different chrom:\s*(\d+)",
    }

    #intiate PE check
    is_paired_end=False

    # Go through files and parse data using regexes
    for f in self.find_log_files('rseqc/bam_stat'):
        d = dict()
        for k, r in regexes.items():
            r_search = re.search(r, f['f'], re.MULTILINE)
            if r_search:
                d[k] = int(r_search.group(1))

        # Calculate some percentages
        if 'total_records' in d:
            t = float(d['total_records'])
            if 'mapq_gte_mapq_cut_unique' in d:
                d['unique_percent'] = (float(d['mapq_gte_mapq_cut_unique']) / t)*100.0
            if 'reads_mapped_in_proper_pairs' in d:
                d['proper_pairs_percent'] = (float(d['reads_mapped_in_proper_pairs']) / t)*100.0

        if len(d) > 0:
            if f['s_name'] in self.bam_stat_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            self.add_data_source(f, section='bam_stat')
            #Check if SE or PE
            if d['read_2'] != 0:
                is_paired_end = True
            self.bam_stat_data[f['s_name']] = d

    # Filter to strip out ignored sample names
    self.bam_stat_data = self.ignore_samples(self.bam_stat_data)

    if len(self.bam_stat_data) > 0:
        # Write to file
        self.write_data_file(self.bam_stat_data, 'multiqc_rseqc_bam_stat')

        # Add to general stats table
        self.general_stats_headers['proper_pairs_percent'] = {
            'title': '% Proper Pairs',
            'description': '% Reads mapped in proper pairs',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn'
        }
        for s_name in self.bam_stat_data:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()

                #Only write if PE, i.e. there is something to write
                if is_paired_end:
                    self.general_stats_data[s_name].update( self.bam_stat_data[s_name] )

        # Make dot plot of counts
        pconfig = {}
        keys = OrderedDict()
        defaults = {
            'min': 0,
            'shared_key': 'read_count',
            'decimalPlaces': 2,
            'modify': lambda x: float(x) / 1000000.0,
        }
        keys['total_records'] = dict(defaults, **{'title': 'Total records' })
        keys['qc_failed'] = dict(defaults, **{'title': 'QC failed' })
        keys['optical_pcr_duplicate'] = dict(defaults, **{'title': 'Duplicates', 'description': 'Optical/PCR duplicate' })
        keys['non_primary_hits'] = dict(defaults, **{'title': 'Non primary hit' })
        keys['unmapped_reads'] = dict(defaults, **{'title': 'Unmapped', 'description': 'Unmapped reads' })
        keys['mapq_lt_mapq_cut_non'] = dict(defaults, **{'title': 'Non-unique', 'description': 'mapq < mapq_cut (non-unique)' })
        keys['mapq_gte_mapq_cut_unique'] = dict(defaults, **{'title': 'Unique', 'description': 'mapq >= mapq_cut (unique)' })
        if is_paired_end:
            keys['read_1'] = dict(defaults, **{'title': 'Read-1' })
            keys['read_2'] = dict(defaults, **{'title': 'Read-2' })
        keys['reads_map_to_sense'] = dict(defaults, **{'title': '+ve strand', 'description': "Reads map to '+'" })
        keys['reads_map_to_antisense'] = dict(defaults, **{'title': '-ve strand', 'description': "Reads map to '-'" })
        keys['non-splice_reads'] = dict(defaults, **{'title': 'Non-splice reads' })
        keys['splice_reads'] = dict(defaults, **{'title': 'Splice reads' })
        if is_paired_end:
            keys['reads_mapped_in_proper_pairs'] = dict(defaults, **{
                        'title': 'Proper pairs', 'description':'Reads mapped in proper pairs' })
            keys['proper-paired_reads_map_to_different_chrom'] = dict(defaults, **{
                        'title': 'Different chrom', 'description': 'Proper-paired reads map to different chrom' })

        self.add_section (
            name = 'Bam Stat',
            anchor = 'rseqc-bam_stat',
            description = 'All numbers reported in millions.',
            plot = beeswarm.plot(self.bam_stat_data, keys)
        )

    # Return number of samples found
    return len(self.bam_stat_data)
