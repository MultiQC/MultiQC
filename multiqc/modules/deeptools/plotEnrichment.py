#!/usr/bin/env python

""" MultiQC submodule to parse output from deepTools plotEnrichment """

import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)

class plotEnrichmentMixin():
    def parse_plotEnrichment(self):
        """Find plotEnrichment output."""
        self.deeptools_plotEnrichment = dict()
        for f in self.find_log_files('deeptools/plotEnrichment'):
            parsed_data = self.parsePlotEnrichment(f['f'], f['fn'])
            for k, v in parsed_data.items():
                if k in self.deeptools_plotEnrichment:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.deeptools_plotEnrichment[k] = v

            if len(parsed_data) > 0:
                self.add_data_source(f, section='plotEnrichment')

        if len(self.deeptools_plotEnrichment) > 0:
            dCounts = {}
            dPercents = {}
            for k, v in self.deeptools_plotEnrichment.items():
                dCounts[k] = dict()
                dPercents[k] = dict()
                for k2, v2 in v.items():
                    dCounts[k][k2] = v2['count']
                    dPercents[k][k2] = v2['percent']
            config = {'data_labels': [
                          {'name': 'Counts', 'ylab': 'Counts in feature'},
                          {'name': 'Percents', 'ylab': 'Percent of reads in feature'}],
                      'ylab': 'Counts',
                      'categories': True,
                      'ymin': 0.0}
            self.add_section(name="plotEnrichment",
                             anchor="plotEnrichment",
                             plot=linegraph.plot([dCounts, dPercents], config))

        return len(self.deeptools_plotEnrichment)

    def parsePlotEnrichment(self, f, fname):
        d = {}
        firstLine = True
        for line in f.splitlines():
            if firstLine:
                firstLine = False
                continue
            cols = line.strip().split("\t")

            if len(cols) != 5:
                log.warning("{} was initially flagged as the output from plotEnrichment, but that seems to not be the case. Skipping...".format(fname))
                return dict()

            if cols[0] not in d:
                d[cols[0]] = dict()
            cols[1] = str(cols[1])
            if cols[1] in d[cols[0]]:
                log.warning("Replacing duplicate sample:featureType {}:{}.".format(cols[0], cols[1]))
            d[cols[0]][cols[1]] = dict()

            try:
                d[cols[0]][cols[1]]["percent"] = float(cols[2])
                d[cols[0]][cols[1]]["count"] = int(cols[3])
            except:
                log.warning("{} was initially flagged as the output from plotEnrichment, but that seems to not be the case. Skipping...".format(fname))
                return dict()
        return d
