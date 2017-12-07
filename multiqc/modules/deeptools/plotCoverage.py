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
            parsed_data = self.parsePlotCoverageStdout(f)
            for k, v in parsed_data.items():
                if k in self.deeptools_plotCoverageStdout:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.deeptools_plotCoverageStdout[k] = v

            if len(parsed_data) > 0:
                self.add_data_source(f, section='plotCoverage')

        self.deeptools_plotCoverageOutRawCounts= dict()
        for f in self.find_log_files('deeptools/plotCoverageOutRawCounts'):
            parsed_data = self.parsePlotCoverageOutRawCounts(f)
            for k, v in parsed_data.items():
                if k in self.deeptools_plotCoverageOutRawCounts:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.deeptools_plotCoverageOutRawCounts[k] = v

            if len(parsed_data) > 0:
                self.add_data_source(f, section='plotCoverage')

        if len(self.deeptools_plotCoverageStdout) > 0:
            header = OrderedDict()
            header["min"] = {
                'title': 'Min',
                'description': 'Minimum Coverage',
                'shared_key': 'coverage'
            }
            header["25%"] = {
                'rid': 'first_quartile',
                'title': '1st Quartile',
                'description': 'First quartile coverage',
                'shared_key': 'coverage'
            }
            header["50%"] = {
                'rid': 'median',
                'title': 'Median',
                'description': 'Median coverage (second quartile)',
                'shared_key': 'coverage'
            }
            header["mean"] = {
                'title': 'Mean',
                'description': 'Mean coverage',
                'shared_key': 'coverage'
            }
            header["75%"] = {
                'rid': 'third_quartile',
                'title': '3rd Quartile',
                'description': 'Third quartile coverage',
                'shared_key': 'coverage'
            }
            header["max"] = {
                'title': 'Max',
                'description': 'Maximum coverage',
                'shared_key': 'coverage'
            }
            header["std"] = {
                'title': 'Std. Dev.',
                'description': 'Coverage standard deviation',
                'shared_key': 'coverage'
            }
            config = {'namespace': 'deepTools plotCoverage'}
            self.add_section(
                name = "Coverage metrics",
                anchor = "deeptools_coverage_metrics",
                plot = table.plot(self.deeptools_plotCoverageStdout, header, config)
            )

        if len(self.deeptools_plotCoverageOutRawCounts) > 0:
            config = {
                'id': 'deeptools_coverage_metrics_plot',
                'title': 'Coverage distribution',
                'xlab': 'Coverage',
                'ylab': 'Fraction of bases sampled'
            }
            self.add_section(
                name = "Coverage distribution",
                anchor = "deeptools_coverage_distribution",
                description = "The fraction of bases with a given number of read/fragment coverage",
                plot = linegraph.plot(self.deeptools_plotCoverageOutRawCounts, config)
            )

        return len(self.deeptools_plotCoverageStdout), len(self.deeptools_plotCoverageOutRawCounts)

    def parsePlotCoverageStdout(self, f):
        d = {}
        firstLine = True
        for line in f['f'].splitlines():
            if firstLine:
                firstLine = False
                continue
            cols = line.strip().split("\t")

            if len(cols) != 8:
                log.warning("{} was initially flagged as the standard output from plotCoverage, but that seems to not be the case. Skipping...".format(f['fn']))
                return dict()

            s_name = self.clean_s_name(cols[0], f['root'])
            if s_name in d:
                log.warning("Replacing duplicate sample {}.".format(s_name))
            d[s_name] = dict()

            try:
                d[s_name]["mean"] = float(cols[1])
                d[s_name]["std"] = float(cols[2])
                d[s_name]["min"] = float(cols[3])
                d[s_name]["25%"] = float(cols[4])
                d[s_name]["50%"] = float(cols[5])
                d[s_name]["75%"] = float(cols[6])
                d[s_name]["max"] = float(cols[7])
            except:
                log.warning("{} was initially flagged as the standard output from plotCoverage, but that seems to not be the case. Skipping...".format(f['fn']))
                return dict()
        return d

    def parsePlotCoverageOutRawCounts(self, f):
        samples = []
        d = {}
        nCols = 0
        nRows = 0
        for line in f['f'].splitlines():
            if line.startswith('#plotCoverage'):
                continue

            cols = line.strip().split('\t')
            if len(cols) < 4:
                log.warning("{} was initially flagged as the output from plotCoverage --outRawCounts, but that seems to not be the case. Skipping...".format(f['fn']))
                return dict()

            if cols[0] == "#'chr'":
                nCols = len(cols)
                for col in cols[3:]:
                    s_name = self.clean_s_name(col.strip("'"), f['root'])
                    samples.append(s_name)
                    d[s_name] = dict()
                continue

            if len(cols) != nCols:
                log.warning("{} was initially flagged as the output from plotCoverage --outRawCounts, but that seems to not be the case. Skipping...".format(f['fn']))
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
