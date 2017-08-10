#!/usr/bin/env python

""" MultiQC module to parse TsTv by alternative allele count from vcftools TsTv-by-count """

import csv
import logging
from collections import defaultdict
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)

class TsTvByCountMixin():

    def parse_tstv_by_count(self):
        """ Create the HTML for the TsTv by alternative allele count linegraph plot. """

        self.vcftools_tstv_by_count = dict()
        for f in self.find_log_files('vcftools/tstv_by_count', filehandles=True):
            d = {}
            for line in f['f'].readlines()[1:]: # don't add the header line (first row)
                key = float(line.split()[0]) # taking the first column (alternative allele count) as key
                val = float(line.split()[3]) # taking Ts/Tv as value
                d[key] = val
            self.vcftools_tstv_by_count[f['s_name']] = d

        # Filter out ignored sample names
        self.vcftools_tstv_by_count = self.ignore_samples(self.vcftools_tstv_by_count)

        if len(self.vcftools_tstv_by_count) == 0:
            return 0

        pconfig = {
            'id': 'vcftools_tstv_by_count_plot',
            'title': 'TsTv by Count',
            'ylab': 'TsTv Ratio',
            'xlab': 'Alternative Allele Count'
        }

        helptext = '''
        `TSTV-BY-COUNT` summarizes the transition to transversion ratio as a function of alternative allele count.
        Note: only using bi-allelic SNPs (see Vcftools's `--TsTv-by-count`).
        '''

        self.add_section(
            name = 'Vcftools TsTv-by-Count',
            anchor = 'vcftools_tstv_by_count',
            description = "Plot of `TSTV-BY-COUNT` - the transition to transversion ratio as a function of alternative allele count from the output of vcftools TsTv-by-count.",
            helptext = helptext,
            plot = linegraph.plot(self.vcftools_tstv_by_count,pconfig)
        )

        return len(self.vcftools_tstv_by_count)

