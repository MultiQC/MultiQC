#!/usr/bin/env python

""" MultiQC submodule to parse output from RSeQC geneBody_coverage.py
http://rseqc.sourceforge.net/#genebody-coverage-py """

from collections import OrderedDict
import logging

from multiqc import config, plots

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
            nrows = 0
            for l in f['f'].splitlines():
                s = l.split()
                if len(keys) == 0:
                    keys = s[1:]
                else:
                    nrows += 1
                    s_name = s[0]
                    if s_name.endswith('.geneBodyCoverage'):
                        s_name = s_name[:-17]
                    if s_name in self.gene_body_cov_hist_counts:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                    self.add_data_source(f, s_name, section='gene_body_coverage')
                    self.gene_body_cov_hist_counts[s_name] = OrderedDict()
                    for k, var in enumerate(s[1:]):
                        self.gene_body_cov_hist_counts[s_name][int(keys[k])] = float(var)
            if nrows == 0:
                log.warning("Empty geneBodyCoverage file found: {}".format(f['fn']))
                
        
        # RSeQC < v2.4
        elif f['f'].startswith('Total reads'):
            if f['s_name'].endswith('.geneBodyCoverage'):
                f['s_name'] = f['s_name'][:-17]
            if f['s_name'] in self.gene_body_cov_hist_counts:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            self.add_data_source(f, section='gene_body_coverage')
            self.gene_body_cov_hist_counts[f['s_name']] = OrderedDict()
            nrows = 0
            for l in f['f'].splitlines():
                s = l.split()
                try:
                    nrows += 1
                    self.gene_body_cov_hist_counts[f['s_name']][int(s[0])] = float(s[1])
                except:
                    pass
            if nrows == 0:
                del self.gene_body_cov_hist_counts[f['s_name']]
                log.warning("Empty geneBodyCoverage file found: {}".format(f['fn']))
    
    if len(self.gene_body_cov_hist_counts) > 0:
        
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
        p_link = '<a href="http://rseqc.sourceforge.net/#genebody-coverage-py" target="_blank">Gene Body Coverage</a>'
        self.sections.append({
            'name': 'Gene Body Coverage',
            'anchor': 'rseqc-gene_body_coverage',
            'content': "<p>"+p_link+" calculates read coverage over gene bodies." \
                " This is used to check if reads coverage is uniform and" \
                " if there is any 5' or 3' bias.</p>" + 
                plots.linegraph.plot([self.gene_body_cov_hist_counts, self.gene_body_cov_hist_percent], pconfig)
        })
    
    # Return number of samples found
    return len(self.gene_body_cov_hist_counts)
    