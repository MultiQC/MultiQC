#!/usr/bin/env python

""" MultiQC module to parse TsTv by alternative allele count from vcftools TsTv-by-count """

import logging
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
            'id': 'vcftools_tstv_by_count',
            'title': 'VCFTools: TsTv by Count',
            'ylab': 'TsTv Ratio',
            'xlab': 'Alternative Allele Count',
            'xmin': 0,
            'ymin': 0,
            'smooth_points': 400, # this limits huge filesizes and prevents browser crashing
            'smooth_points_sumcounts': False
        }

        helptext = '''
        `Transition` is a purine-to-purine or pyrimidine-to-pyrimidine point mutations.
        `Transversion` is a purine-to-pyrimidine or pyrimidine-to-purine point mutation.
        `Alternative allele count` is the number of alternative alleles at the site.
        Note: only bi-allelic SNPs are used (multi-allelic sites and INDELs are skipped.)
        Refer to Vcftools's manual (https://vcftools.github.io/man_latest.html) on `--TsTv-by-count`
        '''

        self.add_section(
            name = 'TsTv by Count',
            anchor = 'vcftools-tstv-by-count',
            description = "Plot of `TSTV-BY-COUNT` - the transition to transversion ratio as a function of alternative allele count from the output of vcftools TsTv-by-count.",
            helptext = helptext,
            plot = linegraph.plot(self.vcftools_tstv_by_count,pconfig)
        )

        return len(self.vcftools_tstv_by_count)

