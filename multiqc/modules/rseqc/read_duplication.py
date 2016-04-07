#!/usr/bin/env python

""" MultiQC submodule to parse output from RSeQC read_duplication.py
http://rseqc.sourceforge.net/#read-duplication-py """

from collections import OrderedDict
import logging
import re

from multiqc import config

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find RSeQC read_duplication reports and parse their data """
    
    # Set up vars
    self.read_dups = dict()
    
    # Go through files and parse data
    for f in self.find_log_files(config.sp['rseqc']['read_duplication_pos']):
        if f['f'].startswith('Occurrence	UniqReadNumber'):
            s_name = f['s_name'].rstrip('.pos.DupRate.xls')
            print('Found {} - "{}"'.format(f['s_name'], s_name))
            self.read_dups[s_name] = OrderedDict()
            for l in f['f'].splitlines():
                s = l.split()
                try:
                    if int(s[0]) <= 500:
                        self.read_dups[s_name][int(s[0])] = int(s[1])
                except:
                    pass
    
    if len(self.read_dups) > 0:
        
        # Log output
        self.sample_count += len(self.read_dups)
        log.info("Found {} read_duplication reports".format(len(self.read_dups)))
        
        # Add line graph to section
        pconfig = {
            'id': 'rseqc_read_dups_plot',
            'title': 'RSeQC: Read Duplication',
            'ylab': 'Number of Reads (log10)',
            'xlab': "Occurance of read",
            'yLog': True,
            'tt_label': "<strong>{point.x} occurances</strong>: {point.y} reads",
        }
        p_link = '<a href="http://rseqc.sourceforge.net/#read-duplication-py" target="_blank">read_duplication.py</a>'
        self.sections.append({
            'name': 'Read Duplication',
            'anchor': 'rseqc-read_dups',
            'content': "<p>"+p_link+" calculates how many alignment positions have a certain number of exact duplicates."\
                " Note - plot truncated at 500 occurances.</p>" + 
                self.plot_xy_data(self.read_dups, pconfig)
        })
    
    
        