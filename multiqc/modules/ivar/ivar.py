#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MultiQC module to parse output files from iVar """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

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
        parsed_data = dict()
        regexes = {
            'total_reads': r'(?:Trimmed primers from )(?:\d+\.\d+\% \()?(\d+)',
            'reads_outside_primer_region': r'^(?:\d+\.\d+\% \()?(\d+)(?:\))?(?:.*[of]?)reads\sstarted',
            'reads_too_short_after_trimming': r'^(?:\d+\.\d+\% \()?(\d+)(?:\))?(?:.*[of]?)reads\swere(?: quality trimmed | shortened)'
        }
        for l in f['f']:
            # Search regexes for stats
            for k, r in regexes.items():
                match = re.search(r, l)
                if match:
                    parsed_data[k] = int(match.group(1))
        return parsed_data


