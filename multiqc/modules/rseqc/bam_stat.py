#!/usr/bin/env python

""" MultiQC submodule to parse output from RSeQC bam_stat.py
http://rseqc.sourceforge.net/#bam-stat-py """

import logging
import re

from multiqc import config

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
    
    # Go through files and parse data using regexes
    for f in self.find_log_files(config.sp['rseqc']['bam_stat']):
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
            self.bam_stat_data[f['s_name']] = d
    
    if len(self.bam_stat_data) > 0:
        
        # Log output
        self.sample_count += len(self.bam_stat_data)
        log.info("Found {} bam_stat reports".format(len(self.bam_stat_data)))
    
        # Write to file
        self.write_data_file(self.bam_stat_data, 'multiqc_rseqc_bam_stat')
        
        # Add to general stats table
        self.general_stats_headers['proper_pairs_percent'] = {
            'title': '% Proper Pairs',
            'description': '% Reads mapped in proper pairs',
            'max': 100,
            'min': 0,
            'scale': 'YlGn',
            'format': '{:.1f}%'
        }
        for s_name in self.bam_stat_data:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update( self.bam_stat_data[s_name] )
        
        # Make dot plot of counts
        # TODO - write dot plot function
    
    
        