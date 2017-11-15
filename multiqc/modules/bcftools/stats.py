#!/usr/bin/env python

""" MultiQC submodule to parse output from Bcftools stats """

import logging
from collections import OrderedDict
from multiqc import config
from multiqc.plots import bargraph, linegraph, table

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
        collapse_complementary = getattr(config, 'bcftools', {}).get('collapse_complementary_changes', False)
        if collapse_complementary:
            types = ['A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T']
        else:
            types = ['A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T',
                     'G>A', 'G>C', 'G>T', 'T>A', 'T>C', 'T>G']

        self.bcftools_stats = dict()
        self.bcftools_stats_indels = dict()
        depth_data = dict()
        for f in self.find_log_files('bcftools/stats'):
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
                    depth_data[s_name] = OrderedDict()
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

                    rc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                    change = s[2].strip()
                    if change not in types:
                        change = '>'.join(rc[n] for n in change.split('>'))

                    field = 'substitution_type_{}'.format(change)
                    value = float(s[3].strip())
                    if field not in self.bcftools_stats[s_name]:
                        self.bcftools_stats[s_name][field] = 0
                    self.bcftools_stats[s_name][field] += value

                # Indel length distributions
                if s[0] == "IDD" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    length = float(s[2].strip())
                    count = float(s[3].strip())
                    self.bcftools_stats_indels[s_name][length] = count

                # Per-sample counts
                if s[0] == "PSC" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    fields = ['variations_hom', 'variations_het']
                    for i, f in enumerate(fields):
                        self.bcftools_stats[s_name][f] = int(s[i + 4].strip())

                # Depth plots
                if s[0] == "DP" and len(s_names) > 0:
                    s_name = s_names[int(s[1])]
                    bin_name = s[2].strip()
                    percent_sites = float(s[-1].strip())
                    depth_data[s_name][bin_name] = percent_sites

        # Filter to strip out ignored sample names
        self.bcftools_stats = self.ignore_samples(self.bcftools_stats)

        if len(self.bcftools_stats) > 0:

            # Write parsed report data to a file
            self.write_data_file(self.bcftools_stats, 'multiqc_bcftools_stats')

            # Stats Table
            stats_headers = self.bcftools_stats_genstats_headers()
            if getattr(config, 'bcftools', {}).get('write_general_stats', True):
                self.general_stats_addcols(self.bcftools_stats, stats_headers, 'Bcftools Stats')
            if getattr(config, 'bcftools', {}).get('write_separate_table', False):
                self.add_section(
                    name='Bcftools Stats',
                    anchor='bcftools-stats',
                    plot=table.plot(self.bcftools_stats, stats_headers))

            # Make bargraph plot of substitution types
            keys = OrderedDict()
            for t in types:
                keys['substitution_type_{}'.format(t)] = {'name': t}
            pconfig = {
                'id': 'bcftools-stats-subtypes',
                'title': 'Bcftools Stats: Substitutions',
                'ylab': '# Substitutions',
                'cpswitch_counts_label': 'Number of Substitutions'
            }
            self.add_section (
                name = 'Variant Substitution Types',
                anchor = 'bcftools-stats',
                plot = bargraph.plot(self.bcftools_stats, keys, pconfig)
            )

            # Make line graph of indel lengths
            if len(self.bcftools_stats_indels) > 0:
                pconfig = {
                    'id': 'bcftools_stats_indel-lengths',
                    'title': 'Bcftools Stats: Indel Distribution',
                    'ylab': 'Count',
                    'xlab': 'InDel Length (bp)',
                    'xDecimals': False,
                    'ymin': 0,
                    # 'tt_label': '<b>{point.x} bp trimmed</b>: {point.y:.0f}',
                }
                self.add_section (
                    name = 'Indel Distribution',
                    anchor = 'bcftools-stats_indel_plot',
                    plot = linegraph.plot(self.bcftools_stats_indels, pconfig)
                )
            # Make line graph of variants per depth
            if len(depth_data) > 0:
                pconfig = {
                    'id': 'bcftools_stats_depth',
                    'title': 'Bcftools Stats: Variant depths',
                    'ylab': 'Fraction of sites (%)',
                    'xlab': 'Variant depth',
                    'ymin': 0,
                    'ymax': 100,
                    'categories': True
                }
                self.add_section (
                    name = 'Variant depths',
                    anchor = 'bcftools-stats_depth_plot',
                    description = 'Read depth support distribution for called variants',
                    plot = linegraph.plot(depth_data, pconfig)
                )

        # Return the number of logs that were found
        return len(self.bcftools_stats)

    def bcftools_stats_genstats_headers(self):
        """ Add key statistics to the General Stats table """
        stats_headers = OrderedDict()
        stats_headers['number_of_records'] = {
            'title': 'Vars',
            'description': 'Variations total',
            'min': 0, 'format': '{:,.0f}',
        }
        stats_headers['variations_hom'] = {
            'title': 'Hom',
            'description': 'Variations homozygous',
            'min': 0, 'format': '{:,.0f}',
        }
        stats_headers['variations_het'] = {
            'title': 'Het',
            'description': 'Variations heterozygous',
            'min': 0, 'format': '{:,.0f}',
        }
        stats_headers['number_of_SNPs'] = {
            'title': 'SNP',
            'description': 'Variation SNPs',
            'min': 0, 'format': '{:,.0f}',
        }
        stats_headers['number_of_indels'] = {
            'title': 'Indel',
            'description': 'Variation Insertions/Deletions',
            'min': 0, 'format': '{:,.0f}',
        }
        stats_headers['tstv'] = {
            'title': 'Ts/Tv',
            'description': 'Variant SNP transition / transversion ratio',
            'min': 0, 'format': '{:,.2f}',
        }
        stats_headers['number_of_MNPs'] = {
            'title': 'MNP',
            'description': 'Variation multinucleotide polymorphisms',
            'min': 0, 'format': '{:,.0f}', "hidden": True,
        }
        return stats_headers
