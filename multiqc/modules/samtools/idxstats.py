#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MultiQC submodule to parse output from Samtools idxstats """

import logging
import re
from collections import OrderedDict
from multiqc import config, plots

# Initialise the logger
log = logging.getLogger(__name__)

class IdxstatsReportMixin():

    def parse_samtools_idxstats(self):
        """ Find Samtools idxstats logs and parse their data """

        self.samtools_idxstats = dict()
        for f in self.find_log_files(config.sp['samtools']['idxstats']):
            parsed_data = parse_single_report(f['f'])
            if len(parsed_data) > 0:
                if f['s_name'] in self.samtools_idxstats:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f, section='idxstats')
                self.samtools_idxstats[f['s_name']] = parsed_data

        if len(self.samtools_idxstats) > 0:

            # Write parsed report data to a file (restructure first)
            self.write_data_file(self.samtools_idxstats, 'multiqc_samtools_idxstats')
            
            # Make a bar plot
            keys = list()
            pdata = dict()
            # Count the total count
            total_mapped = 0
            for s_name in self.samtools_idxstats:
                for k in self.samtools_idxstats[s_name]:
                    total_mapped += self.samtools_idxstats[s_name][k]['mapped']
            for s_name in self.samtools_idxstats:
                pdata[s_name] = OrderedDict()
                for k in self.samtools_idxstats[s_name]:
                    if k not in keys:
                        keys.append(k)
                    # Only add if > 0.1% of total mapped reads
                    mapped = self.samtools_idxstats[s_name][k]['mapped']
                    if mapped >= total_mapped*0.001:
                        pdata[s_name][k] = mapped
            
            pconfig = {
                'id': 'samtools-idxstats-mapped-reads',
                'title': 'Samtools idxstats - Mapped Reads',
                'ylab': '# Reads',
                'categories': keys,
                'tt_label': '<strong>{point.category}:</strong> {point.y}'
            }
            
            self.sections.append({
                'name': 'Samtools idxstats',
                'anchor': 'samtools-idxstats',
                'content': '<p>The <code>samtools idxstats</code> tool counts the number of mapped reads per chromosome. ' +
                    'Chromosomes with &lt; 0.1% of the total aligned reads are omitted from this plot.</p>' +
                    plots.linegraph.plot(pdata, pconfig)
            })
        
        # Return the number of logs that were found
        return len(self.samtools_idxstats)


# idxstats has four columns: chr, length, mapped, unmapped
# http://www.htslib.org/doc/samtools.html

def parse_single_report(f):
    """ Parse a samtools idxstats idxstats """
    
    parsed_data = OrderedDict()
    for l in f.splitlines():
        s = l.split("\t")
        try:
            if int(s[1]) == 0:
                continue
            try:
                pct_mapped = (float(s[2])/(float(s[2])+float(s[3])))*100.0
            except ZeroDivisionError:
                pct_mapped = 0
            parsed_data[s[0]] = {
                'length': int(s[1]),
                'mapped': int(s[2]),
                'unmapped': int(s[3]),
                'pct_mapped': pct_mapped,
                'mapped_per_bp': float(s[2])/float(s[1])
            }
        except (IndexError, ValueError):
            pass
    return parsed_data
    