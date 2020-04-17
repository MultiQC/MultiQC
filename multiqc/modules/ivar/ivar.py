#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MultiQC module to parse output files from iVar """

from __future__ import print_function
from collections import OrderedDict
import logging

#TODO update this
import json
from multiqc.plots import bargraph, linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name = 'ivar',
            anchor = 'ivar',
            href = 'https://github.com/andersen-lab/ivar',
            info = "is a computational package that contains functions broadly useful for viral amplicon-based sequencing."
        )

        # Find and load iVar trim results
        self.ivar_data = dict()
        for f in self.find_log_files('ivar/trim', filehandles=True):
            parsed_data = self.parse_ivar(f)
            if parsed_data is not None and len(parsed_data) > 0:
                self.ivar_data[f['s_name']] = parsed_data
                self.add_data_source(f, f['s_name'])

        # Filter to strip out ignored sample names
        self.ivar_data = self.ignore_samples(self.ivar_data)

        # Warning when no files are found
        if len(self.ivar_data) == 0:
            raise UserWarning

        # Write parsed data to a file
        self.write_data_file(self.ivar_data, 'multiqc_ivar_summary')

        #Found reports or not?
        log.info("Found {} reports".format(len(self.ivar_data)))

        # Report sections TODO
            

    # Parse a ivar report
    def parse_ivar(self, f):
        total_reads = 0
        reads_too_short_after_trimming = 0
        reads_outside_primer_region = 0
        parsed_data = dict()

        for l in f:
            l = l.strip()
            #Find total reads first, iVar 1 vs 2 are different here, thus more complex regex req.
            if "Trimmed primers from" in l:
                match = re.search(r'\(?\d+\)?[^\d+\.\d+\%]', l)
                if match:
                    total_reads = str.replace(match.group(1), '(', '')
                    total_reads = str.replace(total_reads, ')', '')
                    parsed_data['total_reads'] = total_reads
            if "below the minimum length of" in l:
                match = re.search(r'\(?\d{3,}\)?', l)
                if match:
                    reads_too_short_after_trimming = str.replace(match.group(1), '(', '')
                    reads_too_short_after_trimming = str.replace(reads_too_short_after_trimming, ')', '')
                    parsed_data['reads_too_short_after_trimming'] = reads_too_short_after_trimming
            if "reads started outside of primer regions" in l:
                match = re.search(r'\(?\d{3,}\)?', l)
                if match:
                    reads_outside_primer_region = str.replace(match.group(1), '(', '')
                    reads_outside_primer_region = str.replace(reads_outside_primer_region, ')', '')
                    parsed_data['reads_outside_primer_region'] = reads_outside_primer_region
        return parsed_data



