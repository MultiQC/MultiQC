#!/usr/bin/env python

""" MultiQC submodule to parse output from RSeQC junction_annotation.py
http://rseqc.sourceforge.net/#junction-annotation-py """

from collections import OrderedDict
import logging
import re

from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find RSeQC junction_annotation reports and parse their data """

    # Set up vars
    self.junction_annotation_data = dict()
    regexes = {
        'total_splicing_events': r"^Total splicing  Events:\s*(\d+)$",
        'known_splicing_events': r"^Known Splicing Events:\s*(\d+)$",
        'partial_novel_splicing_events': r"^Partial Novel Splicing Events:\s*(\d+)$",
        'novel_splicing_events': r"^Novel Splicing Events:\s*(\d+)$",
        'total_splicing_junctions': r"^Total splicing  Junctions:\s*(\d+)$",
        'known_splicing_junctions': r"^Known Splicing Junctions:\s*(\d+)$",
        'partial_novel_splicing_junctions': r"^Partial Novel Splicing Junctions:\s*(\d+)$",
        'novel_splicing_junctions': r"^Novel Splicing Junctions:\s*(\d+)$",
    }

    # Go through files and parse data using regexes
    for f in self.find_log_files('rseqc/junction_annotation'):
        d = dict()
        for k, r in regexes.items():
            r_search = re.search(r, f['f'], re.MULTILINE)
            if r_search:
                d[k] = int(r_search.group(1))

        # Calculate some percentages
        if 'total_splicing_events' in d:
            t = float(d['total_splicing_events'])
            if 'known_splicing_events' in d:
                d['known_splicing_events_pct'] = (float(d['known_splicing_events']) / t)*100.0
            if 'partial_novel_splicing_events' in d:
                d['partial_novel_splicing_events_pct'] = (float(d['partial_novel_splicing_events']) / t)*100.0
            if 'novel_splicing_events' in d:
                d['novel_splicing_events_pct'] = (float(d['novel_splicing_events']) / t)*100.0
        if 'total_splicing_junctions' in d:
            t = float(d['total_splicing_junctions'])
            if 'known_splicing_junctions' in d:
                d['known_splicing_junctions_pct'] = (float(d['known_splicing_junctions']) / t)*100.0
            if 'partial_novel_splicing_junctions' in d:
                d['partial_novel_splicing_junctions_pct'] = (float(d['partial_novel_splicing_junctions']) / t)*100.0
            if 'novel_splicing_junctions' in d:
                d['novel_splicing_junctions_pct'] = (float(d['novel_splicing_junctions']) / t)*100.0

        if len(d) > 0:
            if f['s_name'] in self.junction_annotation_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            self.add_data_source(f, section='junction_annotation')
            self.junction_annotation_data[f['s_name']] = d

    # Filter to strip out ignored sample names
    self.junction_annotation_data = self.ignore_samples(self.junction_annotation_data)

    if len(self.junction_annotation_data) > 0:

        # Write to file
        self.write_data_file(self.junction_annotation_data, 'multiqc_rseqc_junction_annotation')

        # Plot junction annotations
        keys = [OrderedDict(), OrderedDict()]
        keys[0]['known_splicing_junctions'] = { 'name': 'Known Splicing Junctions' }
        keys[0]['partial_novel_splicing_junctions'] = { 'name': 'Partial Novel Splicing Junctions' }
        keys[0]['novel_splicing_junctions'] = { 'name': 'Novel Splicing Junctions' }
        keys[1]['known_splicing_events'] = { 'name': 'Known Splicing Events' }
        keys[1]['partial_novel_splicing_events'] = { 'name': 'Partial Novel Splicing Events' }
        keys[1]['novel_splicing_events'] = { 'name': 'Novel Splicing Events' }

        pconfig = {
            'id': 'rseqc_junction_annotation_junctions_plot',
            'title': 'RSeQC: Splicing Junctions',
            'ylab': '% Junctions',
            'cpswitch_c_active': False,
            'data_labels': [ 'Junctions', 'Events' ]
        }
        self.add_section (
            name = 'Junction Annotation',
            anchor = 'rseqc_junction_annotation',
            description = '<a href="http://rseqc.sourceforge.net/#junction-annotation-py" target="_blank">Junction annotation</a>' \
                " compares detected splice junctions to" \
                " a reference gene model. An RNA read can be spliced 2" \
                " or more times, each time is called a splicing event.",
            plot = bargraph.plot([self.junction_annotation_data, self.junction_annotation_data], keys, pconfig)
        )

    # Return number of samples found
    return len(self.junction_annotation_data)

