#!/usr/bin/env python

""" MultiQC submodule to parse output from deepTools estimateReadFiltering """

import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.plots import table

# Initialise the logger
log = logging.getLogger(__name__)

class estimateReadFilteringMixin():
    def parse_estimateReadFiltering(self):
        """Find estimateReadFiltering output. Only the output from --table is supported."""
        self.deeptools_estimateReadFiltering = dict()
        for f in self.find_log_files('deeptools/estimateReadFiltering'):
            parsed_data = self.parseEstimateReadFilteringFile(f)
            for k, v in parsed_data.items():
                if k in self.deeptools_estimateReadFiltering:
                    log.warning("Replacing duplicate sample {}.".format(k))
                self.deeptools_estimateReadFiltering[k] = v

            if len(parsed_data) > 0:
                self.add_data_source(f, section='estimateReadFiltering')

        if len(self.deeptools_estimateReadFiltering) > 0:
            header = OrderedDict()
            header["M Entries"] = {
                'title': 'M entries',
                'description': 'Number of entries in the file (millions)'}
            header["pct_Aligned"] = {
                'title': '% Aligned',
                'description': 'Percent of aligned entries',
                'scale': 'YlGn',
                'min': 0,
                'max': 100
            }
            header["pct_Filtered"] = {
                'title': '% Tot. Filtered',
                'description': 'Percent of alignment that would be filtered for any reason.',
                'scale': 'OrRd',
                'min': 0,
                'max': 100
            }
            header["pct_Blacklisted"] = {
                'title': '% Blacklisted',
                'description': 'Percent of alignments falling (at least partially) inside a blacklisted region',
                'scale': 'YlOrRd',
                'min': 0,
                'max': 100
            }
            header["pct_MAPQ"] = {
                'title': '% MAPQ', 'description':
                'Percent of alignments having MAPQ scores below the specified threshold',
                'scale': 'YlOrBn',
                'min': 0,
                'max': 100
            }
            header["pct_Missing_Flags"] = {
                'title': '% Missing Flags',
                'description': 'Percent of alignments lacking at least on flag specified by --samFlagInclude',
                'scale': 'PuRd',
                'min': 0,
                'max': 100
            }
            header["pct_Forbidden_Flags"] = {
                'title': '% Forbidden Flags',
                'description': 'Percent of alignments having at least one flag specified by --samFlagExclude',
                'scale': 'OrRd',
                'min': 0,
                'max': 100
            }
            header["pct_deepTools_Dupes"] = {
                'title': '% deepTools Dupes',
                'description': 'Percent of alignments marked by deepTools as being duplicates',
                'scale': 'PuRd',
                'min': 0,
                'max': 100
            }
            header["pct_Duplication"] = {
                'title': '% Duplication',
                'description': 'Percent of alignments originally marked as being duplicates',
                'scale': 'OrRd',
                'min': 0,
                'max': 100
            }
            header["pct_Singletons"] = {
                'title': '% Singletons',
                 'description': 'Percent of alignments that are singletons (i.e., paired-end reads where the mates don\'t align as a pair',
                 'scale': 'OrRd',
                 'min': 0,
                 'max': 100
                }
            header["pct_Strand_Filtered"] = {
                'title': '% Strand Filtered',
                'description': 'Percent of alignments arising from the wrong strand',
                'scale': 'OrRd',
                'min': 0,
                'max': 100
            }

            tdata = dict()
            for k, v in self.deeptools_estimateReadFiltering.items():
                tdata[k] = {
                    'M Entries': v['total'] / 1000000.0,
                    'pct_Aligned': 100. * v['mapped'] / float(v['total']),
                    'pct_Filtered': 100. * v['filtered'] / float(v['total']),
                    'pct_Blacklisted': 100. * v['blacklisted'] / float(v['total']),
                    'pct_Below_MAPQ': 100. * v['mapq'] / float(v['total']),
                    'pct_Missing_Flags': 100. * v['required flags'] / float(v['total']),
                    'pct_Forbidden_Flags': 100. * v['excluded flags'] / float(v['total']),
                    'pct_deepTools_Dupes': 100. * v['internal dupes'] / float(v['total']),
                    'pct_Duplication': 100. * v['dupes'] / float(v['total']),
                    'pct_Singletons': 100. * v['singletons'] / float(v['total']),
                    'pct_Strand_Filtered': 100. * v['strand'] / float(v['total'])
                }

            config = {'namespace': 'deepTools bamPEFragmentSize'}
            self.add_section(
                name = "Filtering metrics",
                anchor = "estimateReadFiltering",
                description = "Estimated percentages of alignments filtered independently for each setting in `estimateReadFiltering`",
                plot = table.plot(tdata, header, config)
            )

        return len(self.deeptools_estimateReadFiltering)

    def parseEstimateReadFilteringFile(self, f):
        d = {}
        firstLine = True
        for line in f['f'].splitlines():
            if firstLine:
                firstLine = False
                continue
            cols = line.strip().split("\t")

            if len(cols) != 12:
                # This is not really the output from estimateReadFiltering!
                log.warning("{} was initially flagged as the tabular output from estimateReadFiltering, but that seems to not be the case. Skipping...".format(f['fn']))
                return dict()

            s_name = self.clean_s_name(cols[0], f['root'])
            if s_name in d:
                log.debug("Replacing duplicate sample {}.".format(s_name))
            d[s_name] = dict()

            try:
                d[s_name]["total"] = int(cols[1])
                d[s_name]["mapped"] = int(cols[2])
                d[s_name]["blacklisted"] = int(cols[3])
                d[s_name]["filtered"] = float(cols[4])
                d[s_name]["mapq"] = float(cols[5])
                d[s_name]["required flags"] = float(cols[6])
                d[s_name]["excluded flags"] = float(cols[7])
                d[s_name]["internal dupes"] = float(cols[8])
                d[s_name]["dupes"] = float(cols[9])
                d[s_name]["singletons"] = float(cols[10])
                d[s_name]["strand"] = float(cols[11])
            except:
                # Obviously this isn't really the output from estimateReadFiltering
                log.warning("{} was initially flagged as the output from estimateReadFiltering, but that seems to not be the case. Skipping...".format(f['fn']))
                return dict()
        return d
