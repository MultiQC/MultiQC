#!/usr/bin/env python

""" MultiQC module to parse logs from SnpEff """

from __future__ import print_function

import os
from collections import OrderedDict
import logging
import re
from multiqc import config, BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """ SnpEff """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='SnpEff', anchor='snpeff',
                                            href="http://snpeff.sourceforge.net/",
                                            info=" is a tool to annotate SNPs.")

        self.snpeff_data = dict()

        for f in self.find_log_files(config.sp['snpeff'], filehandles=True):
            self.parse_snpeff_log(f)
        
        if len(self.snpeff_data) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning
        
        log.info("Found {} reports".format(len(self.snpeff_data)))
        
        # Write parsed report data to a file
        self.write_data_file(self.snpeff_data, 'multiqc_snpeff')
        
        # General stats table
        self.general_stats()
        
        # Report sections
        self.sections = list()
        self.sections.append({
            'name': 'Counts by effect',
            'anchor': 'snpeff-effects',
            'content': self.count_effects_plot()
        })
    

    def parse_snpeff_log(self, f):
        """ Go through log file looking for snpeff output """
        
        keys = {
            '# Summary table': [
                'Genome', 'Number_of_variants_before_filter', 'Number_of_known_variants',
                'Number_of_effects', 'Genome_total_length', 'Change_rate'
            ],
            '# Effects by impact': [ 'HIGH', 'LOW', 'MODERATE', 'MODIFIER' ],
            '# Effects by functional class': [ 'MISSENSE', 'NONSENSE', 'SILENT', 'Missense_Silent_ratio' ],
            '# Hom/Het table': ['Het', 'Hom', 'Missing'],
            '# Ts/Tv summary': [ 'Transitions', 'Transversions', 'Ts_Tv_ratio' ],
            '# Count by effects': 'all',
            '# Count by genomic region': 'all'
        }
        parsed_data = {}
        section = None
        for l in f['f']:
            l = l.strip()
            if l[:1] == '#':
                section = l
                continue
            s = l.split(',')
            if section in keys:
                if keys[section] == 'all' or any([k in s[0].strip() for k in keys[section]]):
                    try:
                        parsed_data[ s[0].strip() ] = float(s[1].strip())
                    except ValueError:
                        parsed_data[ s[0].strip() ] = s[1].strip()
                    except IndexError:
                        pass
                    if len(s) > 2 and s[2][-1:] == '%':
                        parsed_data[ '{}_percent'.format(s[0].strip()) ] = float(s[2][:-1])
        
        if len(parsed_data) > 0:
            if f['s_name'] in self.snpeff_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            self.add_data_source(f)
            self.snpeff_data[f['s_name']] = parsed_data
    
    def general_stats(self):
        """ Add key SnpEff stats to the general stats table """
        
        headers = OrderedDict()
        headers['Change_rate'] = {
            'title': 'Change rate',
            'scale': 'RdYlBu-rev',
            'min': 0
        }
        headers['Ts_Tv_ratio'] = {
            'title': 'Ts/Tv',
            'description': 'Transitions / Transversions ratio',
        }
        headers['Number_of_variants_before_filter'] = {
            'title': 'M Variants',
            'description': 'Number of variants before filter (millions)',
            'scale': 'PuRd',
            'modify': lambda x: x / 1000000,
            'min': 0
        }
        self.general_stats_addcols(self.snpeff_data, headers)
    
    def count_effects_plot(self):
        """ Generate the HiCUP Truncated reads plot """    
        
        # Choose the categories to ply
        keys = [
            '3_prime_UTR_variant',
            '5_prime_UTR_premature_start_codon_gain_variant',
            '5_prime_UTR_variant',
            'TF_binding_site_variant',
            'downstream_gene_variant',
            'initiator_codon_variant',
            'intergenic_region',
            'intragenic_variant',
            'intron_variant',
            'missense_variant',
            'missense_variant',
            'non_coding_exon_variant',
            'sequence_feature',
            'sequence_feature',
            'splice_region_variant',
            'start_lost',
            'stop_gained',
            'stop_lost',
            'stop_retained_variant',
            'synonymous_variant',
            'upstream_gene_variant',
        ]
        
        # Sort the keys based on the first dataset
        for s_name in self.snpeff_data:
            sorted_keys = sorted(keys, reverse=True, key=self.snpeff_data[s_name].__getitem__)
            break
        
        # Make nicer label names
        pkeys = OrderedDict()
        for k in sorted_keys:
            pkeys[k] = {'name': k.replace('_', ' ').title().replace('Utr', 'UTR') }
        
        # Config for the plot
        pconfig = {
            'title': 'SnpEff: Count by effects',
            'ylab': '# Reads',
        }
        
        return self.plot_bargraph(self.snpeff_data, pkeys, pconfig)
