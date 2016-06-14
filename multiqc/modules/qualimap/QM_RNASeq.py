#!/usr/bin/env python

""" MultiQC Submodule to parse output from Qualimap RNASeq """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config, plots

# Initialise the logger
log = logging.getLogger(__name__)

def parse_reports(self):
    """ Find Qualimap RNASeq reports and parse their data """
    
    sp = config.sp['qualimap']['rnaseq']
    self.qualimap_rnaseq_genome_results = dict()
    regexes = {
        'bam_file': r"bam file\s*=\s*(.+)",
        'reads_aligned': r"read(?:s| pairs) aligned\s*=\s*([\d,]+)",
        'non_unique_alignments': r"non-unique alignments\s*=\s*([\d,]+)",
        'reads_aligned_genes': r"aligned to genes\s*=\s*([\d,]+)",
        'ambiguous_alignments': r"ambiguous alignments\s*=\s*([\d,]+)",
        'not_aligned': r"not aligned\s*=\s*([\d,]+)",
        '5_3_bias': r"5'-3' bias\s*=\s*([\d\.]+)",
        'reads_aligned_exonic': r"exonic\s*=\s*([\d\.]+)",
        'reads_aligned_intronic': r"intronic\s*=\s*([\d\.]+)",
        'reads_aligned_intergenic': r"intergenic\s*=\s*([\d\.]+)",
        'reads_aligned_overlapping_exon': r"overlapping exon\s*=\s*([\d\.]+)",
    }
    for f in self.find_log_files(sp['rnaseq_results']):
        d = dict()
        for k, r in regexes.items():
            r_search = re.search(r, f['f'], re.MULTILINE)
            if r_search:
                try:
                    d[k] = float(r_search.group(1).replace(',',''))
                except ValueError:
                    d[k] = r_search.group(1)
                except UnicodeEncodeError:
                    # Qualimap reports infinity (\u221e) when 3' bias denominator is zero
                    pass
        
        # Check we have an input filename
        if 'bam_file' not in d:
            log.debug("Couldn't find an input filename in genome_results file {}".format(f['fn']))
            return None
        
        # Get a nice sample name
        s_name = self.clean_s_name(d['bam_file'], f['root'])
        
        # Add to general stats table
        for k in ['5_3_bias', 'reads_aligned_genes', 'reads_aligned']:
            try:
                self.general_stats_data[s_name][k] = d[k]
            except KeyError:
                pass
        
        # Save results
        if s_name in self.qualimap_rnaseq_genome_results:
            log.debug("Duplicate genome results sample name found! Overwriting: {}".format(s_name))
        self.qualimap_rnaseq_genome_results[s_name] = d
        self.add_data_source(f, s_name=s_name, section='rna_genome_results')
    
    
    #### Coverage profile
    self.qualimap_rnaseq_cov_hist = dict()
    for f in self.find_log_files(sp['coverage'], filehandles=True):
        s_name = self.get_s_name(f)
        d = dict()
        for l in f['f']:
            if l.startswith('#'):
                continue
            coverage, count = l.split(None, 1)
            coverage = int(round(float(coverage)))
            count = float(count)
            d[coverage] = count
        
        if len(d) == 0:
            log.debug("Couldn't parse contents of coverage histogram file {}".format(f['fn']))
            return None
        
        # Save results
        if s_name in self.qualimap_rnaseq_cov_hist:
            log.debug("Duplicate coverage histogram sample name found! Overwriting: {}".format(s_name))
        self.qualimap_rnaseq_cov_hist[s_name] = d
        self.add_data_source(f, s_name=s_name, section='rna_coverage_histogram')
    
    #### Plots
    
    # Genomic Origin Bar Graph
    if len(self.qualimap_rnaseq_genome_results) > 0:
        gorigin_cats = OrderedDict()
        gorigin_cats['reads_aligned_intergenic'] = {'name': 'Intergenic'}
        gorigin_cats['reads_aligned_intronic'] = {'name': 'Intronic'}
        gorigin_cats['reads_aligned_exonic'] = {'name': 'Exonic'}
        gorigin_pconfig = {
            'title': 'Genomic Origin',
            'cpswitch_c_active': False
        }
        self.sections.append({
            'name': 'Reads genomic origin',
            'anchor': 'qualimap-reads-genomic-origin',
            'content': plots.bargraph.plot(self.qualimap_rnaseq_genome_results, gorigin_cats, gorigin_pconfig)
        })
    
    if len(self.qualimap_rnaseq_cov_hist) > 0:
        self.sections.append({
            'name': 'Coverage Profile Along Genes (total)',
            'anchor': 'qualimap-genome-fraction-coverage',
            'content': plots.linegraph.plot(self.qualimap_rnaseq_cov_hist, {
                'title': 'Coverage Profile Along Genes (total)',
                'ylab': 'Coverage',
                'xlab': 'Transcript Position (bp)',
                'ymin': 0,
                'xmin': 0,
                'xmax': 100,
                'tt_label': '<b>{point.x} bp</b>: {point.y:.0f}%',
            })
        })
    
    
    #### General Stats
    self.general_stats_headers['5_3_bias'] = {
        'title': "5'-3' bias"
    }
    self.general_stats_headers['reads_aligned_genes'] = {
        'title': 'Reads in Genes',
        'description': 'Reads Aligned - Genes (millions)',
        'min': 0,
        'scale': 'PuBu',
        'shared_key': 'read_count',
        'modify': lambda x: x / 1000000,
    }
    self.general_stats_headers['reads_aligned'] = {
        'title': 'Aligned',
        'description': 'Reads Aligned (millions)',
        'min': 0,
        'scale': 'RdBu',
        'shared_key': 'read_count',
        'modify': lambda x: x / 1000000,
    }
    
    # Return the number of reports we found
    return len(self.qualimap_rnaseq_genome_results.keys())
