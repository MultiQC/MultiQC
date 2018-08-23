#!/usr/bin/env python

""" MultiQC submodule to parse output from deepTools plotProfile """

import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)

class plotProfileMixin():
    def parse_plotProfile(self):
        """Find plotProfile output"""
        self.deeptools_plotProfile = dict()
        for f in self.find_log_files('deeptools/plotProfile', filehandles=False):
            parsed_data, bin_labels = self.parsePlotProfileData(f)
            for k, v in parsed_data.items():
                if k in self.deeptools_plotProfile:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.deeptools_plotProfile[k] = v
            if len(parsed_data) > 0:
                self.add_data_source(f, section='plotProfile')

        if len(self.deeptools_plotProfile) > 0:
            config = {
                'id': 'read_distribution_profile',
                'title': 'deeptools: Read Distribution Profile after Annotation',
                'ylab': 'Occurrence',
                'xlab': None,
                'smooth_points': 100,
                'xPlotBands': [
                    {'from': bin_labels.index('TES')+1, 'to': len(bin_labels), 'color': '#f7cfcf'},
                    {'from': bin_labels.index('TSS')+1, 'to': bin_labels.index('TES')+1, 'color': '#ffffe2'},
                    {'from': 1, 'to': bin_labels.index('TSS')+1, 'color': '#e5fce0'},
                ],
                'xPlotLines': [
                    {'width': 1, 'value': bin_labels.index('TES')+1, 'dashStyle': 'Dash', 'color': '#000000'},
                    {'width': 1, 'value': bin_labels.index('TSS')+1, 'dashStyle': 'Dash', 'color': '#000000'},
                ],
            }

            self.add_section (
                name = 'Read Distribution Profile after Annotation',
                anchor = 'read_distribution_profile_plot',
                description="Accumulated view of the distribution of sequence reads related to the closest annotated gene. All annotated genes have been normalized to the same size. Green: {} upstream of gene to {}; Yellow: {} to {}; Pink: {} to {} downstream of gene".format(list(filter(None,bin_labels))[0], list(filter(None,bin_labels))[1], list(filter(None,bin_labels))[1], list(filter(None,bin_labels))[2], list(filter(None,bin_labels))[2], list(filter(None,bin_labels))[3]),
                plot=linegraph.plot(self.deeptools_plotProfile, config)
            )

        return len(self.deeptools_bamPEFragmentSizeDistribution)

    def parsePlotProfileData(self, f):
        d = dict()
        bin_labels = []
        bins = []
        for line in f['f'].splitlines():
            cols = line.rstrip().split("\t")
            if cols[0] == "bin labels":
                for col in cols[2:len(cols)]:
                    if col not in list(filter(None,bin_labels)):
                        bin_labels.append(col)
                    else:
                        break
            elif cols[0] == "bins":
                for col in cols[2:len(cols)]:
                    if len(bins)!=len(bin_labels):
                        bins.append(int(col))
                    else:
                        break
            else:
                s_name = self.clean_s_name(cols[0], f['root'])
                d[s_name] = dict()
                for i in bins:
                    d[s_name].update({i:float(cols[i+1])})

        return d, bin_labels
