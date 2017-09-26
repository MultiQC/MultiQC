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
            parsed_data = self.parseBamPEFile(f)
            for k, v in parsed_data.items():
                if k in self.deeptools_bamPEFragmentSize:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.deeptools_bamPEFragmentSize[k] = v

            if len(parsed_data) > 0:
                self.add_data_source(f, section='bamPEFragmentSize')
        if len(self.deeptools_bamPEFragmentSize) > 0:
            headersSE = OrderedDict()
            headersSE["Reads Sampled"] = {'title': 'Reads Sampled'}
            headersSE["Read Len. Min."] = {'title': 'Read Len. Min.', 'description': 'Minimum read length'}
            headersSE["Read Len. 1st. Qu."] = {'title': 'Read Len. 1st Qu.', 'description': '1st quartile read length'}
            headersSE["Read Len. Mean"] = {'title': 'Read Len. Mean', 'description': 'Mean read length'}
            headersSE["Read Len. Median"] = {'title': 'Read Len. Median', 'description': 'Median read length'}
            headersSE["Read Len. 3rd Qu."] = {'title': 'Read Len. 3rd Qu.', 'description': '3rd quartile read length'}
            headersSE["Read Len. Max"] = {'title': 'Read Len. Max.', 'description': 'Maximum read length'}
            headersSE["Read Len. Std."] = {'title': 'Read Len. Std. Dev.', 'description': 'read length standard deviation'}
            headersSE["Read Med. Abs. Dev."] = {'title': 'Read Len. MAD', 'description': 'read length median absolute deviation'}
            config = {'namespace': 'deepTools bamPEFragmentSize'}
            self.add_section(name="Read lengths",
                             anchor="bamPEFragmentSize",
                             plot=table.plot(self.deeptools_bamPEFragmentSize, headersSE, config))

            headersPE = OrderedDict()
            headersPE["Frag. Sampled"] = {'title': 'Fragments Sampled'}
            headersPE["Frag. Len. Min."] = {'title': 'Fragment Len. Min.', 'description': 'Minimum fragment length'}
            headersPE["Frag. Len. 1st. Qu."] = {'title': 'Fragment Len. 1st Qu.', 'description': '1st quartile fragment length'}
            headersPE["Frag. Len. Mean"] = {'title': 'Fragment Len. Mean', 'description': 'Mean fragment length'}
            headersPE["Frag. Len. Median"] = {'title': 'Fragment Len. Median', 'description': 'Median fragment length'}
            headersPE["Frag. Len. 3rd Qu."] = {'title': 'Fragment Len. 3rd Qu.', 'description': '3rd quartile fragment length'}
            headersPE["Frag. Len. Max"] = {'title': 'Fragment Len. Max.', 'description': 'Maximum fragment length'}
            headersPE["Frag. Len. Std."] = {'title': 'Fragment Len. Std. Dev.', 'description': 'Fragment length standard deviation'}
            headersPE["Frag. Med. Abs. Dev."] = {'title': 'Fragment Len. MAD', 'description': 'Fragment length median absolute deviation'}

            # Are there any PE datasets?
            PE = False
            for k, v in self.deeptools_bamPEFragmentSize.items():
                if 'Frag. Len. Min.' in v:
                    PE = True
                    break
            if PE:
                self.add_section(name="Fragment lengths",
                                 anchor="bamPEFragmentSize",
                                 plot=table.plot(self.deeptools_bamPEFragmentSize, headersPE, config))

            # Read length plot
            config = {'data_labels': [dict(name="Read length distribution", title="Read length distribution", ylab="Read length (bases)"),
                                      dict(name="Fragment length distribution", title="Fragment length distribution", ylab="Fragment length (bases)")],
                      'id': 'bamPEFragmentSize',
                      'title': 'Read/Fragment length distribution',
                      'namespace': 'deepTools bamPEFragmentSize',
                      'ylab': "Read length (bases)",
                      'xlab': "Percentile"}
            SE = dict()
            PE = dict()
            for k, v in self.deeptools_bamPEFragmentSize.items():
                SE[k] = {0: v['Read Len. Min.'],
                         10: v['Read Len. 10%'],
                         20: v['Read Len. 20%'],
                         25: v['Read Len. 1st. Qu.'],
                         30: v['Read Len. 30%'],
                         40: v['Read Len. 40%'],
                         50: v['Read Len. Median'],
                         60: v['Read Len. 60%'],
                         70: v['Read Len. 70%'],
                         75: v['Read Len. 3rd Qu.'],
                         80: v['Read Len. 80%'],
                         90: v['Read Len. 90%'],
                         99: v['Read Len. 99%'],
                         100: v['Read Len. Max']}
                if 'Frag. Len. Min.' not in v:
                    continue
                PE[k] = {0: v['Frag. Len. Min.'],
                         10: v['Frag. Len. 10%'],
                         20: v['Frag. Len. 20%'],
                         25: v['Frag. Len. 1st. Qu.'],
                         30: v['Frag. Len. 30%'],
                         40: v['Frag. Len. 40%'],
                         50: v['Frag. Len. Median'],
                         60: v['Frag. Len. 60%'],
                         70: v['Frag. Len. 70%'],
                         75: v['Frag. Len. 3rd Qu.'],
                         80: v['Frag. Len. 80%'],
                         90: v['Frag. Len. 90%'],
                         99: v['Frag. Len. 99%'],
                         100: v['Frag. Len. Max']}
            self.add_section(name="Read/fragment length distribution",
                             anchor="bamPEFragmentSize",
                             plot=linegraph.plot([SE, PE], config))

        return len(self.deeptools_bamPEFragmentSize)

    def parseBamPEFile(self, f):
        d = {}
        headers = None
        for line in f['f'].splitlines():
            cols = line.rstrip().split("\t")
            if headers is None:
                headers = cols
            else:
                s_name = None
                for idx, h in enumerate(headers):
                    if idx == 0:
                        s_name = self.clean_s_name(cols[0], f['root'])
                        if s_name in d:
                            log.debug("Replacing duplicate sample {}.".format(s_name))
                        d[s_name] = OrderedDict()
                    else:
                        if idx < 19 and cols[1] == "0":
                            # Don't store fragment metrics for SE datasets, they're just 0.
                            continue
                        try:
                            # Most values are ac
                            d[s_name][h] = int(cols[idx])
                        except ValueError:
                            d[s_name][h] = float(cols[idx])

        return d
