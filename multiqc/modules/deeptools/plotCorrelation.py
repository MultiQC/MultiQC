#!/usr/bin/env python

""" MultiQC submodule to parse output from deepTools plotCorrelation """

import logging
import re
from collections import OrderedDict
import numpy as np

from multiqc import config
from multiqc.plots import heatmap

# Initialise the logger
log = logging.getLogger(__name__)

class plotCorrelationMixin():
    def parse_plotCorrelation(self):
        """Find plotCorrelation output"""
        self.deeptools_plotCorrelationData = dict()
        for f in self.find_log_files('deeptools/plotCorrelationData', filehandles=False):
            parsed_data, samples = self.parsePlotCorrelationData(f)
            for k, v in parsed_data.items():
                if k in self.deeptools_plotCorrelationData:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.deeptools_plotCorrelationData[k] = v
            if len(parsed_data) > 0:
                self.add_data_source(f, section='plotCorrelation')

        if len(self.deeptools_plotCorrelationData) > 0:
            config = {
                'id': 'deeptools_correlation_plot',
                'title': 'deeptools: Correlation Plot',
                'square': True,
                'colstops': [
                    [0, '#313695'],
                    [0.1, '#4575b4'],
                    [0.2, '#74add1'],
                    [0.3, '#abd9e9'],
                    [0.4, '#e0f3f8'],
                    [0.5, '#ffffbf'],
                    [0.6, '#fee090'],
                    [0.7, '#fdae61'],
                    [0.8, '#f46d43'],
                    [0.9, '#d73027'],
                    [1, '#a50026'],
                ],
                'decimalPlaces': 2,
                'borderWidth': 0,
                'datalabels': True,
                'datalabel_colour': '<auto>',
            }
            data = []
            for s_name in samples:
                try:
                    data.append(self.deeptools_plotCorrelationData[s_name])
                except KeyError:
                    pass
            if len(data) == 0:
                log.debug('No valid data for correlation plot')
                return None

            self.add_section(
                name="Pairwise Correlation plot",
                anchor="deeptools_correlation",
                description="pairwise correlations",
                plot=heatmap.plot(data, samples, samples, config)
            )

        return len(self.deeptools_plotCorrelationData)

    def parsePlotCorrelationData(self, f):
        d = dict()
        samples = []
        for line in f['f'].splitlines():
            cols = line.split('\t')
            if cols[0] == "#plotCorrelation --outFileCorMatrix":
                continue
            elif cols[0] == "":
                continue
            else:
                c = str(cols[0]).strip("'")
                s_name = self.clean_s_name(c, f['root'])
                samples.append(s_name)
                d[s_name] = []
                for c in cols[1:len(cols)]:
                    d[s_name].append(float(c))
        return d, samples
