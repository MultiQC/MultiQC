#!/usr/bin/env python

""" MultiQC module to parse TsTv by summary output from vcftools TsTv-summary """

import csv
import logging
from collections import defaultdict, OrderedDict
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)

class TsTvSummaryMixin():

    def parse_tstv_summary(self):
        """ Create the HTML for the TsTv summary plot. """

        self.vcftools_tstv_summary = dict()
        for f in self.find_log_files('vcftools/tstv_summary', filehandles=True):
            d = {}
            for line in f['f'].readlines()[1:]: # don't add the header line (first row)
                key = line.split()[0] # taking the first column (MODEL) as key
                val = int(line.split()[1]) # taking the second column (COUNT) as value
                d[key] = val
            self.vcftools_tstv_summary[f['s_name']] = d

        # Filter out ignored sample names
        self.vcftools_tstv_summary = self.ignore_samples(self.vcftools_tstv_summary)

        if len(self.vcftools_tstv_summary) == 0:
            return 0

        # Specifying the categories of the bargraph
        keys = OrderedDict()
        keys['AC'] = { 'name': 'AC' }
        keys['AG'] = { 'name': 'AG' }
        keys['AT'] = { 'name': 'AT' }
        keys['CG'] = { 'name': 'CG' }
        keys['CT'] = { 'name': 'CT' }
        keys['GT'] = { 'name': 'GT' }
        keys['Ts'] = { 'name': 'Ts' }
        keys['Tv'] = { 'name': 'Tv' }

        pconfig = {
            'id': 'vcftools_tstv_summary',
            'title': 'TsTv Summary',
            'ylab': 'Counts',
            'xlab': 'Data Sets'
        }

        helptext = '''
        `TSTV-SUMMARY` summarizes different types of transitions and transversions by count (see Vcftools's `--TsTv-summary`).
        '''

        self.add_section(
            name = 'TsTv Summary',
            anchor = 'vcftools-tstv-summary',
            description = "Plot of `TSTV-SUMMARY` - count of different types of transition and transversion SNPs.",
            helptext = helptext,
            plot = bargraph.plot(self.vcftools_tstv_summary,keys,pconfig)
        )

        return len(self.vcftools_tstv_summary)

