#!/usr/bin/env python

""" MultiQC submodule to parse output from RSeQC geneBody_coverage.py
http://rseqc.sourceforge.net/#genebody-coverage-py """

from collections import OrderedDict
import logging
import re

from multiqc import config

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find RSeQC gene_body_coverage reports and parse their data """
    
    # Set up vars
    self.gene_body_cov_hist_counts = dict()
    self.gene_body_cov_hist_percent = dict()
    
    # TODO - Do separate parsing step to find skewness values
    # and add these to the general stats table?
    
    # Go through files and parse data
    for f in self.find_log_files(config.sp['rseqc']['gene_body_coverage']):
        
        # RSeQC >= v2.4
        if f['f'].startswith('Percentile'):
            keys = []
            for l in f['f'].splitlines():
                s = l.split()
                if len(keys) == 0:
                    keys = s[1:]
                else:
                    s_name = s[0]
                    self.gene_body_cov_hist_counts[s_name] = OrderedDict()
                    for k, var in enumerate(s[1:]):
                        self.gene_body_cov_hist_counts[s_name][int(keys[k])] = float(var)
        
        # RSeQC < v2.4
        elif f['f'].startswith('Total reads'):
            self.gene_body_cov_hist_counts[f['s_name']] = OrderedDict()
            for l in f['f'].splitlines():
                s = l.split()
                try:
                    self.gene_body_cov_hist_counts[f['s_name']][int(s[0])] = float(s[1])
                except:
                    pass
    
    if len(self.gene_body_cov_hist_counts) > 0:
        
        # Log output
        self.sample_count += len(self.gene_body_cov_hist_counts)
        log.info("Found {} geneBody_coverage reports".format(len(self.gene_body_cov_hist_counts)))
        
        # Make a normalised percentage version of the data
        for s_name in self.gene_body_cov_hist_counts:
            self.gene_body_cov_hist_percent[s_name] = OrderedDict()
            total = sum( self.gene_body_cov_hist_counts[s_name].values() )
            for k, v in self.gene_body_cov_hist_counts[s_name].items():
                self.gene_body_cov_hist_percent[s_name][k] = (v/total)*100
        
        # Add line graph to section
        pconfig = {
            'id': 'rseqc_gene_body_coverage_plot',
            'title': 'RSeQC: Gene Body Coverage',
            'ylab': 'Coverage',
            'xlab': "Gene Body Percentile (5' -> 3')",
            'xmin': 0,
            'xmax': 100,
            'tt_label': "<strong>{point.x}% from 5'</strong>: {point.y:.2f}",
            'data_labels': [
                {'name': 'Counts', 'ylab': 'Coverage'},
                {'name': 'Percentages', 'ylab': 'Percentage Coverage'}
            ]
        }
        self.sections.append({
            'name': 'Gene Body Coverage',
            'anchor': 'rseqc-gene_body_coverage',
            'content': "<p>Read coverage over gene body. This module" \
                " is used to check if reads coverage is uniform and" \
                " if there is any 5' or 3' bias.</p>" + 
                self.plot_xy_data([self.gene_body_cov_hist_counts, self.gene_body_cov_hist_percent], pconfig)
        })
    
    
        