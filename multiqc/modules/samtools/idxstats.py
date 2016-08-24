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
            
            # Prep the data for the plots
            keys = list()
            pdata = dict()
            xy_counts = dict()
            # Count the total count
            total_mapped = 0
            for s_name in self.samtools_idxstats:
                for k in self.samtools_idxstats[s_name]:
                    total_mapped += self.samtools_idxstats[s_name][k]['mapped']
            for s_name in self.samtools_idxstats:
                pdata[s_name] = OrderedDict()
                xy_counts[s_name] = dict()
                x_count = False
                y_count = False
                # Loop through each chromosome
                for k in self.samtools_idxstats[s_name]:
                    if k not in keys:
                        keys.append(k)
                    # Only add if > 0.1% of total mapped reads
                    mapped = self.samtools_idxstats[s_name][k]['mapped']
                    if mapped >= total_mapped*0.001:
                        pdata[s_name][k] = mapped
                    # Save counts if chromosome x or y
                    if k.lower() == 'x' or k.lower() == 'chrx':
                        x_count = mapped
                    if k.lower() == 'y' or k.lower() == 'chry':
                        y_count = mapped
                # Only save these counts if we have both x and y
                if x_count and y_count:
                    xy_counts[s_name]['x'] = x_count
                    xy_counts[s_name]['y'] = y_count
            
            # X/Y ratio plot
            if len(xy_counts) > 0:
                xy_keys = OrderedDict()
                xy_keys['x'] = { 'name': 'Chromosome X' }
                xy_keys['y'] = { 'name': 'Chromosome Y' }
                pconfig = {
                    'id': 'samtools-idxstats-xy-plot',
                    'title': 'Samtools idxstats - chrXY mapped reads',
                    'ylab': 'Percent of X+Y Reads',
                    'cpswitch_counts_label': 'Number of Reads',
                    'cpswitch_percent_label': 'Percent of X+Y Reads',
                    'cpswitch_c_active': False
                }
                self.sections.append({
                    'name': 'XY counts',
                    'anchor': 'samtools-idxstats-xy-counts',
                    'content': plots.bargraph.plot(xy_counts, xy_keys, pconfig)
                })
                
            
            # Mapped reads per chr line plot
            pconfig = {
                'id': 'samtools-idxstats-mapped-reads-plot',
                'title': 'Samtools idxstats - Mapped reads per contig',
                'ylab': '# Reads',
                'xlab': 'Chromosome Name',
                'categories': keys,
                'tt_label': '<strong>{point.category}:</strong> {point.y}'
            }
            
            self.sections.append({
                'name': 'Samtools idxstats',
                'anchor': 'samtools-idxstats',
                'content': '<p>The <code>samtools idxstats</code> tool counts the number of mapped reads per chromosome / contig. ' +
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
    