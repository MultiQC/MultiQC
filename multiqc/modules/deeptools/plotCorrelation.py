#!/usr/bin/env python

""" MultiQC submodule to parse output from deepTools plotCorrelation """

import logging

from multiqc.plots import heatmap

# Initialise the logger
log = logging.getLogger(__name__)


class plotCorrelationMixin:
    def parse_plotCorrelation(self):
        """Find plotCorrelation output"""
        self.deeptools_plotCorrelationData = dict()
        for f in self.find_log_files("deeptools/plotCorrelationData", filehandles=False):
            parsed_data, samples = self.parsePlotCorrelationData(f)
            for k, v in parsed_data.items():
                if k in self.deeptools_plotCorrelationData:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.deeptools_plotCorrelationData[k] = v
            if len(parsed_data) > 0:
                self.add_data_source(f, section="plotCorrelation")

        self.deeptools_plotCorrelationData = self.ignore_samples(self.deeptools_plotCorrelationData)

        if len(self.deeptools_plotCorrelationData) > 0:
            config = {
                "id": "deeptools_correlation_plot",
                "title": "deeptools: Correlation Plot",
            }
            data = []
            for s_name in samples:
                try:
                    data.append(self.deeptools_plotCorrelationData[s_name])
                except KeyError:
                    pass
            if len(data) == 0:
                log.debug("No valid data for correlation plot")
                return None

            self.add_section(
                name="Correlation heatmap",
                anchor="deeptools_correlation",
                description="Pairwise correlations of samples based on distribution of sequence reads",
                plot=heatmap.plot(data, samples, samples, config),
            )

        return len(self.deeptools_plotCorrelationData)

    def parsePlotCorrelationData(self, f):
        d = dict()
        samples = []
        for line in f["f"].splitlines():
            cols = line.split("\t")
            if cols[0] == "#plotCorrelation --outFileCorMatrix":
                continue
            elif cols[0] == "":
                continue
            else:
                c = str(cols[0]).strip("'")
                s_name = self.clean_s_name(c, f)
                samples.append(s_name)
                d[s_name] = []
                for c in cols[1 : len(cols)]:
                    d[s_name].append(float(c))
        return d, samples
