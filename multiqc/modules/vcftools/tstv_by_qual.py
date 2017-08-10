#!/usr/bin/env python

""" MultiQC module to parse TsTv by quality output from vcftools TsTv-by-qual """

import csv
import logging
from collections import defaultdict
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)

class TsTvByQualMixin():

    def parse_tstv_by_qual(self):
        """ Create the HTML for the TsTv by quality linegraph plot. """

        self.vcftools_tstv_by_qual = dict()
        for f in self.find_log_files('vcftools/tstv_by_qual', filehandles=True):
            d = {}
            for line in f['f'].readlines()[1:]: # don't add the header line (first row)
                key = float(line.split()[0]) # taking the first column (QUAL_THRESHOLD) as key
                val = float(line.split()[6]) # taking Ts/Tv_GT_QUAL_THRESHOLD as value
                d[key] = val
            self.vcftools_tstv_by_qual[f['s_name']] = d

        # Filter out ignored sample names
        self.vcftools_tstv_by_qual = self.ignore_samples(self.vcftools_tstv_by_qual)

        if len(self.vcftools_tstv_by_qual) == 0:
            return 0

        pconfig = {
            'id': 'vcftools_tstv_by_qual_plot',
            'title': 'TsTv by Qual',
            'ylab': 'TsTv Ratio',
            'xlab': 'SNP Quality Threshold'
        }

        helptext = '''
        `TSTV-BY-QUAL` summarizes the transition to transversion ratio as a function of SNP quality threshold.
        Note: only using bi-allelic SNPs (see Vcftools's `--TsTv-by-qual`).
        '''

        self.add_section(
            name = 'Vcftools TsTv-by-Qual',
            anchor = 'vcftools_tstv_by_qual',
            description = "Plot of `TSTV-BY-QUAL` - the transition to transversion ratio as a function of SNP quality from the output of vcftools TsTv-by-qual.",
            helptext = helptext,
            plot = linegraph.plot(self.vcftools_tstv_by_qual,pconfig)
        )

        return len(self.vcftools_tstv_by_qual)

