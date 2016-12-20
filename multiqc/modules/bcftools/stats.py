#!/usr/bin/env python

""" MultiQC submodule to parse output from Bcftools stats """

import logging
from collections import OrderedDict
from multiqc import config, plots

# Initialise the logger
log = logging.getLogger(__name__)


class StatsReportMixin():
    """ Mixin loaded by the bcftools MultiqcModule class """

    def parse_bcftools_stats(self):
        """
        Find bcftools stats logs and parse their data
          Bcftools stats reports contain 'sets' of data, which can
          have multiple vcf files each (but usually don't). Here,
          we treat each 'set' as a MultiQC sample, taking the first
          input filename for each set as the name.
        """

        self.bcftools_stats = dict()
        self.bcftools_stats_indels = dict()
        for f in self.find_log_files(config.sp['bcftools']['stats']):
            s_names = list()
            for line in f['f'].splitlines():
                s = line.split("\t")
                # Get the sample names - one per 'set'
                if s[0] == "ID":
                    s_name = self.clean_s_name(s[2], f['root'])
                    s_names.append(s_name)
                    if s_name in self.bcftools_stats:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                    self.add_data_source(f, s_name, section='stats')
                    self.bcftools_stats[s_name] = dict()
                    self.bcftools_stats_indels[s_name] = dict()
                    self.bcftools_stats_indels[s_name][0] = None # Avoid joining line across missing 0

                # Parse key stats
                if s[0] == "SN" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    field = s[2].strip()[:-1]
                    field = field.replace(' ', '_')
                    value = float(s[3].strip())
                    self.bcftools_stats[s_name][field] = value

                # Parse transitions/transversions stats
                if s[0] == "TSTV" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    fields = ['ts', 'tv', 'tstv', 'ts_1st_ALT', 'tv_1st_ALT', 'tstv_1st_ALT']
                    for i, f in enumerate(fields):
                        value = float(s[i+2].strip())
                        self.bcftools_stats[s_name][f] = value

                # Parse substitution types
                if s[0] == "ST" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    field = 'substitution_type_{}'.format(s[2].strip())
                    value = float(s[3].strip())
                    self.bcftools_stats[s_name][field] = value

                # Indel length distributions
                if s[0] == "IDD" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    length = float(s[2].strip())
                    count = float(s[3].strip())
                    self.bcftools_stats_indels[s_name][length] = count

        if len(self.bcftools_stats) > 0:

            # Write parsed report data to a file
            self.write_data_file(self.bcftools_stats, 'multiqc_bcftools_stats')

            # General Stats Table
            self.bcftools_stats_genstats_table()

            # Make bargraph plot of substitution types
            types = ['A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G']
            keys = OrderedDict()
            for t in types:
                keys['substitution_type_{}'.format(t)] = {'name': t}
            pconfig = {
                'id': 'bcftools-stats-subtypes',
                'title': 'Bcftools Stats: Substitutions',
                'ylab': '# Substitutions',
                'cpswitch_counts_label': 'Number of Substitutions'
            }
            self.sections.append({
                'name': 'Variant Substitution Types',
                'anchor': 'bcftools-stats',
                'content': plots.bargraph.plot(self.bcftools_stats, keys, pconfig)
            })

            # Make line graph of indel lengths
            if len(self.bcftools_stats_indels) > 1:
                pconfig = {
                    'id': 'bcftools_stats_indel-lengths',
                    'title': 'Bcftools Stats: Indel Distribution',
                    'ylab': 'Count',
                    'xlab': 'InDel Length (bp)',
                    'xDecimals': False,
                    'ymin': 0,
                    # 'tt_label': '<b>{point.x} bp trimmed</b>: {point.y:.0f}',
                }
                self.sections.append({
                    'name': 'Indel Distribution',
                    'anchor': 'bcftools-stats',
                    'content': plots.linegraph.plot(self.bcftools_stats_indels, pconfig)
                })

        # Return the number of logs that were found
        return len(self.bcftools_stats)

    def bcftools_stats_genstats_table(self):
        """ Add key statistics to the General Stats table """
        stats_headers = OrderedDict()
        stats_headers['number_of_SNPs'] = {
            'title': 'M SNPs',
            'description': 'Number of SNPs (millions)',
            'min': 0,
            'format': '{:.2f}',
            'modify': lambda x: x / 1000000
        }
        stats_headers['tstv'] = {
            'title': 'ts/tv',
            'description': 'SNP transitions / transversions ratio',
            'min': 0,
            'format': '{:.2f}',
        }
        stats_headers['number_of_indels'] = {
            'title': 'M Indels',
            'description': 'Number of indels (millions)',
            'min': 0,
            'format': '{:.2f}',
            'modify': lambda x: x / 1000000
        }
        stats_headers['number_of_MNPs'] = {
            'title': 'M MNPs',
            'description': 'Number of MNPs (millions)',
            'min': 0,
            'format': '{:.2f}',
            'modify': lambda x: x / 1000000
        }
        stats_headers['number_of_others'] = {
            'title': 'M Others',
            'description': 'Number of others (millions)',
            'min': 0,
            'format': '{:.2f}',
            'modify': lambda x: x / 1000000
        }
        self.general_stats_addcols(self.bcftools_stats, stats_headers, 'Bcftools Stats')
