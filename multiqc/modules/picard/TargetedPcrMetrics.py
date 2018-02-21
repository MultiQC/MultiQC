#!/usr/bin/env python

""" MultiQC submodule to parse output from Picard TargetedPcrMetrics """

from collections import OrderedDict
import logging
import os
import re

from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """ Find Picard TargetedPcrMetrics reports and parse their data """

    # Set up vars
    self.picard_pcrmetrics_data = dict()
    self.picard_pcrmetrics_samplestats = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files('picard/pcr_metrics', filehandles=True):
        s_name = None
        for l in f['f']:
            # New log starting
            if 'TargetedPcrMetrics' in l and 'INPUT' in l:
                s_name = None
                # Pull sample name from input
                fn_search = re.search(r"INPUT(?:=|\s+)(\[?[^\s]+\]?)", l, flags=re.IGNORECASE)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1).strip('[]'))
                    s_name = self.clean_s_name(s_name, f['root'])

            if s_name is not None:
                if 'TargetedPcrMetrics' in l and '## METRICS CLASS' in l:
                    keys = f['f'].readline().strip("\n").split("\t")
                    vals = f['f'].readline().strip("\n").split("\t")
                    if len(vals) == len(keys):
                        if s_name in self.picard_pcrmetrics_data:
                            log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
                        self.add_data_source(f, s_name, section='TargetedPcrMetrics')
                        self.picard_pcrmetrics_data[s_name] = dict()
                        for i, k in enumerate(keys):
                            try:
                                # Multiply percentages by 100
                                if k.startswith('PCT_'):
                                    vals[i] = float(vals[i]) * 100.0
                                self.picard_pcrmetrics_data[s_name][k] = float(vals[i])
                            except ValueError:
                                self.picard_pcrmetrics_data[s_name][k] = vals[i]

    # Filter to strip out ignored sample names
    self.picard_pcrmetrics_data = self.ignore_samples(self.picard_pcrmetrics_data)

    if len(self.picard_pcrmetrics_data) > 0:

        # Write parsed data to a file
        self.write_data_file(self.picard_pcrmetrics_data, 'multiqc_picard_pcrmetrics')

        # Add to general stats table
        self.general_stats_headers['PCT_AMPLIFIED_BASES'] = {
            'title': '% Amplified Bases',
            'description': 'The fraction of aligned bases that mapped to or near an amplicon.',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'scale': 'BrBG'
        }
        self.general_stats_headers['MEDIAN_TARGET_COVERAGE'] = {
            'title': 'Median Target Coverage',
            'description': 'The median coverage of reads that mapped to target regions of an experiment.',
            'min': 0,
            'suffix': 'X',
            'scale': 'GnBu',
        }

        for s_name in self.picard_pcrmetrics_data:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update( self.picard_pcrmetrics_data[s_name] )

        # Bar plot of ignored bases
        keys = OrderedDict()
        keys['ON_AMPLICON_BASES'] = {'name': 'On-amplicon bases'}
        keys['NEAR_AMPLICON_BASES'] = {'name': 'Near-amplicon bases'}
        keys['OFF_AMPLICON_BASES'] = {'name': 'Off-amplicon bases', 'color': '#f28f43'}

        # Config for the plot
        pconfig = {
            'id': 'picard_pcr_metrics_bases',
            'title': 'Picard: PCR Amplicon Bases',
            'ylab': '# Bases',
            'cpswitch_counts_label': '# Bases',
            'hide_zero_cats': False
        }

        self.add_section (
            name = 'PCR Amplicon Bases',
            anchor = 'picard-pcrmetrics-bases',
            description = 'Metrics about reads obtained from targeted PCR experiments.',
            helptext = '''
            This plot shows the number of bases aligned on or near to amplified regions of the genome.

            * `ON_AMPLICON_BASES`: The number of `PF_BASES_ALIGNED` that mapped to an amplified region of the genome.
            * `NEAR_AMPLICON_BASES`: The number of `PF_BASES_ALIGNED` that mapped to within a fixed interval of an amplified region, but not on a baited region.
            * `OFF_AMPLICON_BASES`: The number of `PF_BASES_ALIGNED` that mapped neither on or near an amplicon.

            For more information see the [Picard documentation](https://broadinstitute.github.io/picard/picard-metric-definitions.html#TargetedPcrMetrics).''',
            plot = bargraph.plot(self.picard_pcrmetrics_data, keys, pconfig)
        )


    # Return the number of detected samples to the parent module
    return len(self.picard_pcrmetrics_data)

