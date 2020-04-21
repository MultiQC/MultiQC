#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MultiQC module to parse output files from iVar """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import bargraph, heatmap
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name = 'iVar',
            anchor = 'iVar',
            href = 'https://github.com/andersen-lab/ivar',
            info = "is a computational package that contains functions broadly useful for viral amplicon-based sequencing."
        )

        # Find and load iVar trim results
        self.ivar_data = dict()
        self.parsed_primers = dict()
        for f in self.find_log_files('ivar/trim', filehandles=True):
            parsed_data = self.parse_ivar(f)
            if parsed_data is not None and len(parsed_data) > 0:
                self.ivar_data[f['s_name']] = parsed_data
                self.add_data_source(f, section='trimming')
        for f in self.find_log_files('ivar/trim'):
            parsed_primers = self.parse_ivar_primer_stats(f)
            if parsed_primers is not None and len(parsed_primers) > 0:
                self.parsed_primers[f['s_name']] = parsed_primers

        # Filter to strip out ignored sample names
        self.ivar_data = self.ignore_samples(self.ivar_data)

        # Warning when no files are found
        if len(self.ivar_data) == 0:
            raise UserWarning

        # Write parsed data to a file
        self.write_data_file(self.ivar_data, 'multiqc_ivar_summary')
        
        #Primers too
        self.write_data_file(self.parsed_primers, 'multiqc_ivar_primers')
        
        #Found reports or not?
        log.info("Found {} reports".format(len(self.ivar_data)))

        # Basic Stats Table
        self.ivar_general_stats_table()

        #Heatmap info
        self.status_heatmap()
        
    # Parse an ivar report
    def parse_ivar(self, f):
        parsed_data = dict()
        regexes = {
            'mapped_reads': r'(?:Found\s)(\d+)(?:\smapped)',
            'total_reads': r'(?:Trimmed primers from )(?:\d+\.\d+\% \()?(\d+)',
            'reads_outside_primer_region': r'^(?:\d+\.\d+\% \()?(\d+)(?:\))?(?:.*[of]?)reads\s(?:that\s)?started',
            'reads_too_short_after_trimming': r'^(?:\d+\.\d+\% \()?(\d+)(?:\))?(?:.*[of]?)reads\swere(?: quality trimmed | shortened)'
        }
        for l in f['f']:
            # Search regexes for stats
            for k, r in regexes.items():
                match = re.search(r, l)
                if match:
                    parsed_data[k] = int(match.group(1))
        
        return parsed_data

    #Parse Primer stats appropriately
    def parse_ivar_primer_stats(self, f):
        primers = OrderedDict()
        regex = "^(.*)(?:\t+)(\d+$)"
        # Search regexes for stats
        for l in f['f'].splitlines():
            match = re.search(regex, l)
            if match:
                primer = match.group(1)
                counts = int(match.group(2))
                primers[primer] = counts
        log.info("Found {} primers".format(len(primers)))
        return primers

    # Add to general stats table

    def ivar_general_stats_table(self):
        """ Take the parsed stats from the iVAR report and add it to the
        basic stats table"""

        headers = OrderedDict()
        headers['mapped_reads'] = {
            'title': 'Total mapped reads',
            'description': 'Total number of mapped reads in iVar input.',
            'min': 0,
            'scale': 'RdYlGn-rev',
            'format': '{:,.0f}'
        }
        headers['total_reads'] = {
            'title': 'Primer trimmed reads',
            'description': 'Total number of reads where primer trimming was performed.',
            'min': 0,
            'scale': 'PuRd',
            'format': '{:,.0f}'
        }
        headers['reads_too_short_after_trimming'] = {
            'title': 'Fail minimum length reads',
            'description': 'Number of reads too short (<30bp) after primer trimming',
            'min': 0,
            'scale': 'RdYlGn',
            'format': '{:,.0f}'
        }
        headers['reads_outside_primer_region'] = {
            'title': 'Reads outside primer region',
            'description': 'Number of reads outside the primer region',
            'scale': 'RdYlGn-rev',
            'format': '{:,.0f}'
        }
        self.general_stats_addcols(self.ivar_data, headers)
    

    def status_heatmap(self):
        """ Heatmap showing information on each primer found for every sample """
        #Top level dict contains sample IDs + OrderedDict(primer, counts)
        
        final_data = list()
        final_xcats = list()
        final_ycats = list()

        for k,v in self.parsed_primers.items():
            final_ycats.append(k)
            tmp_prim_val = list()
            for prim,val in v.items():
                final_xcats.add(prim)
                tmp_prim_val.append(val)
            final_data.append(tmp_prim_val)

        if data is not None:
            pconfig = {
                'id': 'ivar-primer-count-heatmap',
                'decimalPlaces': 0,
                'square': False, 
                'title': 'iVar: Number of primers found for each sample'
            }

            self.add_section (
                name = 'iVar Primer Counts',
                plot = heatmap.plot(final_data, final_xcats, final_ycats, pconfig)
            )
