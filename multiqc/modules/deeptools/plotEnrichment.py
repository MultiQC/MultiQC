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
            parsed_data = self.parsePlotEnrichment(f)
            for k, v in parsed_data.items():
                if k in self.deeptools_plotEnrichment:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.deeptools_plotEnrichment[k] = v

            if len(parsed_data) > 0:
                self.add_data_source(f, section='plotEnrichment')

        if len(self.deeptools_plotEnrichment) > 0:
            dCounts = OrderedDict()
            dPercents = OrderedDict()
            for sample, v in self.deeptools_plotEnrichment.items():
                dCounts[sample] = OrderedDict()
                dPercents[sample] = OrderedDict()
                for category, v2 in v.items():
                    dCounts[sample][category] = v2['count']
                    dPercents[sample][category] = v2['percent']
            config = {'data_labels': [
                          {'name': 'Counts in features', 'ylab': 'Counts in feature'},
                          {'name': 'Percents in features', 'ylab': 'Percent of reads in feature'}],
                      'id': 'deeptools_enrichment_plot',
                      'title': 'Signal enrichment per feature',
                      'ylab': 'Counts in feature',
                      'categories': True,
                      'ymin': 0.0}
            self.add_section(name="Signal enrichment per feature",
                             description="Signal enrichment per feature according to plotEnrichment",
                             anchor="deeptools_enrichment",
                             plot=linegraph.plot([dCounts, dPercents], pconfig=config))

        return len(self.deeptools_plotEnrichment)

    def parsePlotEnrichment(self, f):
        d = {}
        firstLine = True
        for line in f['f'].splitlines():
            if firstLine:
                firstLine = False
                continue
            cols = line.strip().split("\t")

            if len(cols) != 5:
                log.warning("{} was initially flagged as the output from plotEnrichment, but that seems to not be the case. Skipping...".format(f['fn']))
                return dict()

            s_name = self.clean_s_name(cols[0], f['root'])
            if s_name not in d:
                d[s_name] = dict()
            cols[1] = str(cols[1])
            if cols[1] in d[s_name]:
                log.warning("Replacing duplicate sample:featureType {}:{}.".format(s_name, cols[1]))
            d[s_name][cols[1]] = dict()

            try:
                d[s_name][cols[1]]["percent"] = float(cols[2])
                d[s_name][cols[1]]["count"] = int(cols[3])
            except:
                log.warning("{} was initially flagged as the output from plotEnrichment, but that seems to not be the case. Skipping...".format(f['fn']))
                return dict()
        return d
