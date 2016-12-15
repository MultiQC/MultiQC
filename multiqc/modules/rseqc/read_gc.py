#!/usr/bin/env python

""" MultiQC submodule to parse output from RSeQC read_GC.py
http://rseqc.sourceforge.net/#read-gc-py """

from collections import OrderedDict
import logging

from multiqc import config, plots, plots

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find RSeQC read_GC reports and parse their data """

    # Set up vars
    self.read_gc = dict()
    self.read_gc_pct = dict()

    # Go through files and parse data
    for f in self.find_log_files(config.sp['rseqc']['read_gc']):

        if f['f'].startswith('GC%	read_count'):
            gc = list()
            counts = list()
            for l in f['f'].splitlines():
                s = l.split()
                try:
                    gc.append(float(s[0]))
                    counts.append(float(s[1]))
                except:
                    pass
            if len(gc) > 0:
                sorted_gc_keys = sorted(range(len(gc)), key=lambda k: gc[k])
                total = sum(counts)
                if f['s_name'] in self.read_gc:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f, section='read_GC')
                self.read_gc[f['s_name']] = OrderedDict()
                self.read_gc_pct[f['s_name']] = OrderedDict()
                for i in sorted_gc_keys:
                    self.read_gc[f['s_name']][gc[i]] = counts[i]
                    self.read_gc_pct[f['s_name']][gc[i]] = (counts[i]/total)*100

    if len(self.read_gc) > 0:

        # Add line graph to section
        pconfig = {
            'id': 'rseqc_read_gc_plot',
            'title': 'RSeQC: Gene Body Coverage',
            'ylab': 'Number of Reads',
            'xlab': "GC content (%)",
            'xmin': 0,
            'xmax': 100,
            'tt_label': "<strong>{point.x}% GC</strong>: {point.y:.2f}",
            'data_labels': [
                {'name': 'Counts', 'ylab': 'Number of Reads'},
                {'name': 'Percentages', 'ylab': 'Percentage of Reads'}
            ]
        }
        p_link = '<a href="http://rseqc.sourceforge.net/#read-gc-py" target="_blank">read_GC</a>'
        self.sections.append({
            'name': 'Read GC Content',
            'anchor': 'rseqc-read_gc',
            'content': "<p>"+p_link+" calculates a histogram of read GC content.</p>" +
                plots.linegraph.plot([self.read_gc, self.read_gc_pct], pconfig)
        })

    # Return number of samples found
    return len(self.read_gc)



