#!/usr/bin/env python

""" MultiQC submodule to parse output from RSeQC infer_experiment.py 
http://rseqc.sourceforge.net/#infer-experiment-py """

from collections import OrderedDict
import logging
import re

from multiqc import config, plots

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find RSeQC infer_experiment reports and parse their data """
    
    # Set up vars 
    self.infer_exp = dict()
    regexes = {
        'pe_sense': r"\"1\+\+,1--,2\+-,2-\+\": (\d\.\d*)",
        'pe_antisense': r"\"1\+-,1-\+,2\+\+,2--\": (\d\.\d*)",
        'se_sense': r"\"\+\+,--\": (\d\.\d*)",
        'se_antisense': r"\+-,-\+\": (\d\.\d*)",
        'failed': r"\"Fraction of reads failed to determine: (\d\.\d*)"
    }
    
    # Go through files and parse data using regexes
    for f in self.find_log_files(config.sp['rseqc']['infer_experiment']):
        d = dict()
        for k, r in regexes.items():
            r_search = re.search(r, f['f'], re.MULTILINE)
            if r_search:
                d[k] = float(r_search.group(1))
                
        if len(d) > 0:
            if f['s_name'] in self.infer_exp:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            self.add_data_source(f, section='infer_experiment')
            self.infer_exp[f['s_name']] = d
    
    if len(self.infer_exp) > 0:
        
        # Write to file
        self.write_data_file(self.infer_exp, 'multiqc_rseqc_infer_experiment')
        
        # Plot bar graph of groups
        keys = OrderedDict()
        keys['pe_sense'] = {'name': "Paired End Sense"}
        keys['pe_antisense'] = {'name': "Paired End Antisense"}
        keys['se_sense'] = {'name': "Single End Sense"}
        keys['se_antisense'] = {'name': "Single End Antisense"}
        keys['failed'] = {'name': "Undetermined"}
        # Config for the plot
        pconfig = {
            'id': 'rseqc_infer_experiment_plot',
            'title': 'RSeQC: Infer experiment',
            'ylab': '# Tags',
            'cpswitch': False,
            'cpswitch_counts_label': 'Number of Tags',
        }
        
        p_link = '<a href="http://rseqc.sourceforge.net/#infer-experiment-py" target="_blank">Infer experiment</a>'
        self.sections.append({
            'name': 'Infer experiment',
            'anchor': 'rseqc-infer_experiment',
            'content': "<p>"+p_link+" counts how reads and read pairs match the strandedness of overlapping transcripts. It can be used to infer whether RNA-seq library preps are stranded (sense or antisense) .</p>" + 
                plots.bargraph.plot(self.infer_exp, keys, pconfig)
        })
    # Return number of samples found
    return len(self.infer_exp)
    
