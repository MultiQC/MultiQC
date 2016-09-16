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
    first_regexes = {
        'Paired_single': r"This is ([A-Za-z]*) Data"
    }
    second_regexes = {
        'pe_sense': r"\"1\+\+,1--,2\+-,2-\+\": (\d\.\d*)",
        'pe_antisense': r"\"1\+-,1-\+,2\+\+,2--\": (\d\.\d*)",
        'se_sense': r"\"\+\+,--\": (\d\.\d*)",
        'se_antisense': r"\+-,-\+\": (\d\.\d*)",
        'failed': r"\"Fraction of reads failed to determine: (\d\.\d*)"
    }
    
    # Go through files and parse data using regexes
    for f in self.find_log_files(config.sp['rseqc']['infer_experiment']):
        d = dict()
        for k, r in first_regexes.items():
            r_search = re.search(r, f['f'], re.MULTILINE)
            if r_search:
                d[k] = str(r_search.group(1))

        for k, r in second_regexes.items():
            r_search = re.search(r, f['f'], re.MULTILINE)
            if r_search:
                d[k] = float(r_search.group(1))
                
               # d['{}_total_bases'.format(k)] = int(r_search.group(1))
               # d['{}_tag_count'.format(k)] = int(r_search.group(2))
               # d['{}_tags_kb'.format(k)] = float(r_search.group(2))
        
        # Calculate some percentages for parsed file
        #if 'total_tags' in d:
        #    t = float(d['total_tags'])
        #    pcts = dict()
        #    for k in d:
        #        if k.endswith('_tag_count'):
        #            pk = '{}_tag_pct'.format(k[:-10])
        #            pcts[pk] = (float(d[k]) / t)*100.0
        #    d.update(pcts)
        
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
        keys['pe_sense'] = {'name': "PairedEnd Sense"}
        keys['pe_antisense'] = {'name': "PairedEnd Antisense"}
        keys['se_sense'] = {'name': "SingleEnd Sense"}
        keys['se_antisense'] = {'name': "SingleEnd antisense"}
        keys['failed'] = {'name': "Undetermined"}
        # Config for the plot   ----  redo
        pconfig = {
            'id': 'rseqc_infer_experiment_plot',
            'title': 'RSeQC: Read Distribution',
            'ylab': '# Tags',
            'cpswitch': False,
            'cpswitch_counts_label': 'Number of Tags',
            'cpswitch_c_active': False
        }
        
        p_link = '<a href="http://rseqc.sourceforge.net/#infer-experiment-py" target="_blank">Infer experiment</a>'
        self.sections.append({
            'name': 'Infer experiment',
            'anchor': 'rseqc-infer_experiment',
            'content': "<p>"+p_link+" Infer experiment is used to 'guess' how RNA-seq sequencing were configured,in particulary how reads were stranded for strand-specific RNA-seq data, through comparing the 'strandness of reads with the 'strandness of transcripts'.</p>" + 
                plots.bargraph.plot(self.infer_exp, keys, pconfig)
        })
    print self.infer_exp 
    # Return number of samples found
    return len(self.infer_exp)
    
