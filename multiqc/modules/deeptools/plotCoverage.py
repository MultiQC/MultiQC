#!/usr/bin/env python

""" MultiQC submodule to parse output from deepTools plotCoverage """

import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.plots import table, linegraph

# Initialise the logger
log = logging.getLogger(__name__)

class plotCoverageMixin():
    def parse_plotCoverage(self):
        """Find plotCoverage output. Both stdout and --outRawCounts"""
        self.deeptools_plotCoverageStdout = dict()
        for f in self.find_log_files('deeptools/plotCoverageStdout'):
            parsed_data = self.parsePlotCoverageStdout(f['f'], f['fn'])
            for k, v in parsed_data.items():
                if k in self.deeptools_plotCoverageStdout:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.deeptools_plotCoverageStdout[k] = v

            if len(parsed_data) > 0:
                self.add_data_source(f, section='plotCoverage')

        self.deeptools_plotCoverageOutRawCounts= dict()
        for f in self.find_log_files('deeptools/plotCoverageOutRawCounts'):
            parsed_data = self.parsePlotCoverageOutRawCounts(f['f'], f['fn'])
            for k, v in parsed_data.items():
                if k in self.deeptools_plotCoverageOutRawCounts:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.deeptools_plotCoverageOutRawCounts[k] = v

            if len(parsed_data) > 0:
                self.add_data_source(f, section='plotCoverage')

        if len(self.deeptools_plotCoverageStdout) > 0:
            header = OrderedDict()
            header["mean"] = {'title': 'Mean Cov.', 'description': 'Mean coverage'}
            header["std"] = {'title': 'Std. Dev.', 'description': 'Coverage standard deviation'}
            header["min"] = {'title': 'Min. Cov.', 'description': 'Minimum Coverage'}
            header["25%"] = {'title': '1st Q.', 'description': 'First quartile coverage'}
            header["50%"] = {'title': 'Med. Cov.', 'description': 'Median coverage (second quartile)'}
            header["75%"] = {'title': '3rd Q.', 'description': 'Third quartile coverage'}
            header["max"] = {'title': 'Max. Cov.', 'description': 'Maximum coverage'}
            self.add_section(name="plotCoverage (standard output)",
                             anchor="plotCoverage",
                             plot=table.plot(self.deeptools_plotCoverageStdout, header))

        if len(self.deeptools_plotCoverageOutRawCounts) > 0:
            config = dict(xlab='Coverage', ylab='Fraction of bases sampled')
            self.add_section(name="plotCoverage",
                             anchor="plotCoverage",
                             plot=linegraph.plot(self.deeptools_plotCoverageOutRawCounts, config))



        return len(self.deeptools_plotCoverageStdout), len(self.deeptools_plotCoverageOutRawCounts)

    def parsePlotCoverageStdout(self, f, fname):
        d = {}
        firstLine = True
        for line in f.splitlines():
            if firstLine:
                firstLine = False
                continue
            cols = line.strip().split("\t")

            if len(cols) != 8:
                log.warning("{} was initially flagged as the standard output from plotCoverage, but that seems to not be the case. Skipping...".format(fname))
                return dict()

            if cols[0] in d:
                log.warning("Replacing duplicate sample {}.".format(cols[0]))
            d[cols[0]] = dict()

            try:
                d[cols[0]]["mean"] = float(cols[1])
                d[cols[0]]["std"] = float(cols[2])
                d[cols[0]]["min"] = float(cols[3])
                d[cols[0]]["25%"] = float(cols[4])
                d[cols[0]]["50%"] = float(cols[5])
                d[cols[0]]["75%"] = float(cols[6])
                d[cols[0]]["max"] = float(cols[7])
            except:
                log.warning("{} was initially flagged as the standard output from plotCoverage, but that seems to not be the case. Skipping...".format(fname))
                return dict()
        return d

    def parsePlotCoverageOutRawCounts(self, f, fname):
        samples = []
        d = {}
        nCols = 0
        nRows = 0
        for line in f.splitlines():
            cols = line.strip().split('\t')
            if len(cols) < 4:
                log.warning("{} was initially flagged as the output from plotCoverage --outRawCounts, but that seems to not be the case. Skipping...".format(fname))
                return dict()

            if cols[0] == "#'chr'":
                nCols = len(cols)
                for col in cols[3:]:
                    samples.append(col.strip("'"))
                    d[col.strip("'")] = dict()
                continue

            if len(cols) != nCols:
                log.warning("{} was initially flagged as the output from plotCoverage --outRawCounts, but that seems to not be the case. Skipping...".format(fname))
                return dict()

            for i, v in enumerate(cols[3:]):
                v = float(v)
                if v not in d[samples[i]]:
                    d[samples[i]][v] = 0
                d[samples[i]][v] += 1
            nRows += 1

        # Convert values to a fraction
        nRows = float(nRows)
        for k, v in d.items():
            for k2, v2 in v.items():
                v[k2] /= nRows

        return d
