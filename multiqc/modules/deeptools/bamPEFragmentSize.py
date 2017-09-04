#!/usr/bin/env python

""" MultiQC submodule to parse output from deepTools bamPEFragmentSize """

import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.plots import table, linegraph

# Initialise the logger
log = logging.getLogger(__name__)

class bamPEFragmentSizeMixin():
    def parse_bamPEFragmentSize(self):
        """Find bamPEFragmentSize output. Only the output from --table is supported."""
        self.deeptools_bamPEFragmentSize = dict()
        for f in self.find_log_files('deeptools/bamPEFragmentSize'):
            parsed_data = self.parseBamPEFile(f['f'], f['fn'])
            for k, v in parsed_data.items():
                if k in self.deeptools_bamPEFragmentSize:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.deeptools_bamPEFragmentSize[k] = v

            if len(parsed_data) > 0:
                self.add_data_source(f, section='bamPEFragmentSize')
        if len(self.deeptools_bamPEFragmentSize) > 0:
            headersSE = OrderedDict()
            headersSE["Reads Sampled"] = {'title': 'Reads Sampled'}
            headersSE["Read Length Minimum"] = {'title': 'Read Len. Min.', 'description': 'Minimum read length'}
            headersSE["Read Length 1st Quartile"] = {'title': 'Read Len. 1st Qu.', 'description': '1st quartile read length'}
            headersSE["Read Length Mean"] = {'title': 'Read Len. Mean', 'description': 'Mean read length'}
            headersSE["Read Length Median"] = {'title': 'Read Len. Median', 'description': 'Median read length'}
            headersSE["Read Length 3rd Quartile"] = {'title': 'Read Len. 3rd Qu.', 'description': '3rd quartile read length'}
            headersSE["Read Length Maximum"] = {'title': 'Read Len. Max.', 'description': 'Maximum read length'}
            headersSE["Read Length Std. Dev."] = {'title': 'Read Len. Std. Dev.', 'description': 'read length standard deviation'}
            headersSE["Read Length Med. Abs. Dev."] = {'title': 'Read Len. MAD', 'description': 'read length median absolute deviation'}
            self.add_section(name="bamPEFragmentSize (read length metrics)",
                             anchor="bamPEFragmentSize",
                             plot=table.plot(self.deeptools_bamPEFragmentSize, headersSE))

            headersPE = OrderedDict()
            headersPE["Fragments Sampled"] = {'title': 'Fragments Sampled'}
            headersPE["Fragment Length Minimum"] = {'title': 'Fragment Len. Min.', 'description': 'Minimum read length'}
            headersPE["Fragment Length 1st Quartile"] = {'title': 'Fragment Len. 1st Qu.', 'description': '1st quartile read length'}
            headersPE["Fragment Length Mean"] = {'title': 'Fragment Len. Mean', 'description': 'Mean read length'}
            headersPE["Fragment Length Median"] = {'title': 'Fragment Len. Median', 'description': 'Median read length'}
            headersPE["Fragment Length 3rd Quartile"] = {'title': 'Fragment Len. 3rd Qu.', 'description': '3rd quartile read length'}
            headersPE["Fragment Length Maximum"] = {'title': 'Fragment Len. Max.', 'description': 'Maximum read length'}
            headersPE["Fragment Length Std. Dev."] = {'title': 'Fragment Len. Std. Dev.', 'description': 'read length standard deviation'}
            headersPE["Fragment Length Med. Abs. Dev."] = {'title': 'Fragment Len. MAD', 'description': 'read length median absolute deviation'}

            # Are there any PE datasets?
            PE = False
            for k, v in self.deeptools_bamPEFragmentSize.items():
                if 'Fragment Length Minimum' in v:
                    PE = True
                    break
            if PE:
                self.add_section(name="bamPEFragmentSize (fragment length metrics)",
                                 anchor="bamPEFragmentSize",
                                 plot=table.plot(self.deeptools_bamPEFragmentSize, headersPE))

            # Read length plot
            config = {'data_labels': [dict(name="Read length distribution", title="Read length distribution", ylab="Read length (bases)"),
                                      dict(name="Fragment length distribution", title="Fragment length distribution", ylab="Fragment length (bases)")],
                      'ylab': "Read length (bases)",
                      'xlab': "Percentile"}
            SE = dict()
            PE = dict()
            for k, v in self.deeptools_bamPEFragmentSize.items():
                SE[k] = {0: v['Read Length Minimum'],
                         10: v['R10'],
                         20: v['R20'],
                         25: v['Read Length 1st Quartile'],
                         30: v['R30'],
                         40: v['R40'],
                         50: v['Read Length Median'],
                         60: v['R60'],
                         70: v['R70'],
                         75: v['Read Length 3rd Quartile'],
                         80: v['R80'],
                         90: v['R90'],
                         99: v['R99'],
                         100: v['Read Length Maximum']}
                if 'Fragment Length Minimum' not in v:
                    continue
                PE[k] = {0: v['Fragment Length Minimum'],
                         10: v['F10'],
                         20: v['F20'],
                         25: v['Fragment Length 1st Quartile'],
                         30: v['F30'],
                         40: v['F40'],
                         50: v['Fragment Length Median'],
                         60: v['F60'],
                         70: v['F70'],
                         75: v['Fragment Length 3rd Quartile'],
                         80: v['F80'],
                         90: v['F90'],
                         99: v['F99'],
                         100: v['Fragment Length Maximum']}
            self.add_section(name="bamPEFragmentSize (read/fragment length distribution)",
                             anchor="bamPEFragmentSize",
                             plot=linegraph.plot([SE, PE], config))

        return len(self.deeptools_bamPEFragmentSize)

    def parseBamPEFile(self, f, fname):
        d = {}
        firstLine = True
        for line in f.splitlines():
            if firstLine:
                firstLine = False
                continue
            cols = line.strip().split("\t")

            if len(cols) != 37:
                # This is not really the output from bamPEFragmentSize!
                log.warning("{} was initially flagged as the tabular output from bamPEFragmentSize, but that seems to not be the case. Skipping...".format(fname))
                return dict()

            if cols[0] in d:
                log.warning("Replacing duplicate sample {}.".format(cols[0]))
            d[cols[0]] = dict()

            try:
                d[cols[0]]["Reads Sampled"] = int(cols[19])
                d[cols[0]]["Read Length Minimum"] = float(cols[20])
                d[cols[0]]["Read Length 1st Quartile"] = float(cols[21])
                d[cols[0]]["Read Length Mean"] = float(cols[22])
                d[cols[0]]["Read Length Median"] = float(cols[23])
                d[cols[0]]["Read Length 3rd Quartile"] = float(cols[24])
                d[cols[0]]["Read Length Maximum"] = float(cols[25])
                d[cols[0]]["Read Length Std. Dev."] = float(cols[26])
                d[cols[0]]["Read Length Med. Abs. Dev."] = float(cols[27])

                # Not in the table
                d[cols[0]]["R10"] = float(cols[28])
                d[cols[0]]["R20"] = float(cols[29])
                d[cols[0]]["R30"] = float(cols[30])
                d[cols[0]]["R40"] = float(cols[31])
                d[cols[0]]["R60"] = float(cols[32])
                d[cols[0]]["R70"] = float(cols[33])
                d[cols[0]]["R80"] = float(cols[34])
                d[cols[0]]["R90"] = float(cols[35])
                d[cols[0]]["R99"] = float(cols[36])
                if cols[1] != "0":
                    d[cols[0]]["Fragments Sampled"] = int(cols[1])
                    d[cols[0]]["Fragment Length Minimum"] = float(cols[2])
                    d[cols[0]]["Fragment Length 1st Quartile"] = float(cols[3])
                    d[cols[0]]["Fragment Length Mean"] = float(cols[4])
                    d[cols[0]]["Fragment Length Median"] = float(cols[5])
                    d[cols[0]]["Fragment Length 3rd Quartile"] = float(cols[6])
                    d[cols[0]]["Fragment Length Maximum"] = float(cols[7])
                    d[cols[0]]["Fragment Length Std. Dev."] = float(cols[8])
                    d[cols[0]]["Fragment Length Med. Abs. Dev."] = float(cols[9])

                    # Not in the table
                    d[cols[0]]["F10"] = float(cols[10])
                    d[cols[0]]["F20"] = float(cols[11])
                    d[cols[0]]["F30"] = float(cols[12])
                    d[cols[0]]["F40"] = float(cols[13])
                    d[cols[0]]["F60"] = float(cols[14])
                    d[cols[0]]["F70"] = float(cols[15])
                    d[cols[0]]["F80"] = float(cols[16])
                    d[cols[0]]["F90"] = float(cols[17])
                    d[cols[0]]["F99"] = float(cols[18])
            except:
                # Obviously this isn't really the output from bamPEFragmentSize
                log.warning("{} was initially flagged as the tabular output from bamPEFragmentSize, but that seems to not be the case. Skipping...".format(fname))
                return dict()
        return d
