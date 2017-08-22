#!/usr/bin/env python

""" MultiQC submodule to parse output from RSeQC infer_experiment.py
http://rseqc.sourceforge.net/#infer-experiment-py """

from collections import OrderedDict
import logging
import re

from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find RSeQC infer_experiment reports and parse their data """

    # Set up vars
    self.infer_exp = dict()
    regexes = {
        'pe_sense': r"\"1\+\+,1--,2\+-,2-\+\": (\d\.\d+)",
        'pe_antisense': r"\"1\+-,1-\+,2\+\+,2--\": (\d\.\d+)",
        'se_sense': r"\"\+\+,--\": (\d\.\d+)",
        'se_antisense': r"\+-,-\+\": (\d\.\d+)",
        'failed': r"Fraction of reads failed to determine: (\d\.\d+)"
    }

    # Go through files and parse data using regexes
    for f in self.find_log_files('rseqc/infer_experiment'):
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

    # Filter to strip out ignored sample names
    self.infer_exp = self.ignore_samples(self.infer_exp)

    if len(self.infer_exp) > 0:

        # Write to file
        self.write_data_file(self.infer_exp, 'multiqc_rseqc_infer_experiment')

        # Merge PE and SE for plot
        pdata = dict()
        for s_name, vals in self.infer_exp.items():
            pdata[s_name] = dict()
            for k, v in vals.items():
                v *= 100.0 # Multiply to get percentage
                if k[:2] == 'pe' or k[:2] == 'se':
                    k = k[3:]
                pdata[s_name][k] = v + pdata[s_name].get(k, 0)

        # Plot bar graph of groups
        keys = OrderedDict()
        keys['sense'] = {'name': "Sense"}
        keys['antisense'] = {'name': "Antisense"}
        keys['failed'] = {'name': "Undetermined"}

        # Config for the plot
        pconfig = {
            'id': 'rseqc_infer_experiment_plot',
            'title': 'RSeQC: Infer experiment',
            'ylab': '% Tags',
            'ymin': 0,
            'ymax': 100,
            'tt_percentages': False,
            'ylab_format': '{value}%',
            'cpswitch': False
        }

        self.add_section (
            name = 'Infer experiment',
            anchor = 'rseqc-infer_experiment',
            description = '<a href="http://rseqc.sourceforge.net/#infer-experiment-py" target="_blank">Infer experiment</a>' \
                " counts the percentage of reads and read pairs that match the strandedness of overlapping transcripts." \
                " It can be used to infer whether RNA-seq library preps are stranded (sense or antisense).",
            plot = bargraph.plot(pdata, keys, pconfig)
        )

    # Return number of samples found
    return len(self.infer_exp)

