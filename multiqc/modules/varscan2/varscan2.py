#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MultiQC module to parse output files from VarScan2 """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name = 'VarScan2',
            anchor = 'VarScan2',
            href = 'http://dkoboldt.github.io/varscan/',
            info = "variant detection in massively parallel sequencing data"
        )

        # Find and load VarScan2 reports - there are 3 different ones, but all with identical content (differentiated by header)
        self.varscan2_data = dict()
        for f in self.find_log_files('varscan2/mpileup2snp', filehandles=True):
            parsed_data = self.parse_varscan(f)
            s_name = self.clean_s_name(parsed_data['sample_name'], f['root'])
            #Drop existing sample_name from dict now
            del parsed_data['sample_name']
            if parsed_data is not None and len(parsed_data) > 0:
                self.varscan2_data[s_name] = parsed_data
                self.add_data_source(f, s_name)
        for f in self.find_log_files('varscan2/mpileup2indel', filehandles=True):
            parsed_data = self.parse_varscan(f)
            s_name = self.clean_s_name(parsed_data['sample_name'], f['root'])
            #Drop existing sample_name from dict now
            del parsed_data['sample_name']
            if parsed_data is not None and len(parsed_data) > 0:
                self.varscan2_data[s_name] = parsed_data
                self.add_data_source(f, s_name)
        for f in self.find_log_files('varscan2/mpileup2cns', filehandles=True):
            parsed_data = self.parse_varscan(f)
            s_name = self.clean_s_name(parsed_data['sample_name'], f['root'])
            #Drop existing sample_name from dict now
            del parsed_data['sample_name']
            if parsed_data is not None and len(parsed_data) > 0:
                self.varscan2_data[s_name] = parsed_data
                self.add_data_source(f, s_name)

        # Filter to strip out ignored sample names
        self.varscan2_data = self.ignore_samples(self.varscan2_data)

        # Warning when no files are found
        if len(self.varscan2_data) == 0:
            raise UserWarning

        # Write parsed data to a file
        self.write_data_file(self.varscan2_data, 'multiqc_varscan2_summary')

        #Found reports or not?
        log.info("Found {} reports".format(len(self.varscan2_data)))

        # Basic Stats Table
        self.varscan2_general_stats_table()

        #Basic barplot section
        self.add_section(
            description = 'This plot shows the total number of variant positions detected, grouped by variant type.',
            plot = self.varscan2_report_plot()
        )
        #Reported variants plot
        self.add_section(
            description = 'This plot shows the total number of *reported* variant positions, grouped by variant type.',
            plot = self.varscan2_reported_plot()
        )

    # Varscan2 reports in SNP mode only snps, in indel only indel and in CNS all variants found
    #Total variants = SNPs + Indels
    # Parse a VarScan2 report
    def parse_varscan(self, f):
        parsed_data = dict()
        regexes = {
            'sample_name': r'(?:Reading input from )(\w+.+)',
            'min_coverage': r'(?:Min coverage:)\s(\d+)',
            'min_reads2': r'(?:Min reads2:)\s(\d+)',
            'min_var_freq': r'(?:Min var freq:)\s(\d+\.\d+)',
            'min_avg_qual': r'(?:Min avg qual:)\s(\d+)',
            'p_val_threshold': r'(?:P-value thresh:)\s(\d+\.\d+)',
            'bases_in_pileup': r'(\d+)(?:\sbases in pileup file)',
            'variant_pos_total': r'(\d+)(?:\svariant positions \()',
            'variant_pos_snps': r'(?:\d+)(?:\svariant positions \()(\d+)',
            'variant_pos_indels': r'(?:\d+)(?:\svariant positions \()(?:\d+ SNP, )(\d+)',
            'bases_failed_strand_filter': r'(\d+)(?:\swere failed by the strand-filter)',
            'variant_pos_reported_total': r'(\d+)(?:\svariant positions reported)',
            'variant_pos_reported_snps': r'(?:\d+)(?:\svariant positions reported \()(\d+)',
            'variant_pos_reported_indels': r'(?:\d+)(?:\svariant positions reported \()(?:\d+ SNP, )(\d+)'
        }
        for l in f['f']:
            # Search regexes for stats
            for k, r in regexes.items():
                match = re.search(r, l)
                if match:
                    if k not in ['sample_name','p_val_threshold','min_var_freq']:
                        parsed_data[k] = int(match.group(1))
                    if k == 'sample_name':
                        parsed_data[k] = match.group(1)
                    if k in ['p_val_threshold','min_var_freq']:
                        parsed_data[k] = float(match.group(1))
        return parsed_data

    # Add to general stats table
    def varscan2_general_stats_table(self):
        """ Take the parsed stats from the VarScan2 report and add it to the
        basic stats table"""

        headers = OrderedDict()
        headers['min_coverage'] = {
            'title': 'Minimum coverage',
            'description': 'Minimum coverage at variant position observed',
            'min': 0,
            'scale': 'RdYlGn-rev',
            'hidden': True,
            'format': '{:,.0f}'
        }
        headers['min_reads2'] = {
            'title': 'Minimum reads2',
            'description': 'Whatever this means',
            'min': 0,
            'scale': 'RdYlBu',
            'hidden': True,
            'format': '{:,.0f}'
        }
        headers['min_var_freq'] = {
            'title': 'Minium variant frequency',
            'description': 'Minimum variant frequency of variant observed',
            'min': 0,
            'hidden': True,
            'scale': 'Set2',
            'format': '{:,.0f}'
        }
        headers['min_avg_qual'] = {
            'title': 'Minium average quality',
            'description': 'Minimum average quality of variant observed',
            'min': 0,
            'hidden': True,
            'scale': 'BuGn',
            'format': '{:,.0f}'
        }
        headers['p_val_threshold'] = {
            'title': 'P-Value Threshold',
            'description': 'P-Value Threshold for keeping variants',
            'min': 0,
            'hidden': True,
            'scale': 'Set2',
            'format': '{:,.0f}'
        }
        headers['bases_in_pileup'] = {
            'title': 'Bases in Pileup',
            'description': 'Number of bases in pileup input for VarScan2',
            'min': 0,
            'hidden': True,
            'scale': 'Greens',
            'format': '{:,.0f}'
        }
        headers['variant_pos_total'] = {
            'title': 'Total variants detected',
            'description': 'Total number of variants detected',
            'min': 0,
            'scale': 'BrBg',
            'format': '{:,.0f}',
            'shared_key': 'total_variants'
        }
        headers['variant_pos_snps'] = {
            'title': 'Total SNPs detected',
            'description': 'Total number of SNPs detected',
            'min': 0,
            'scale': 'BrBg',
            'format': '{:,.0f}',
            'shared_key': 'total_variants'
        }
        headers['variant_pos_indels'] = {
            'title': 'Total INDELs detected',
            'description': 'Total number of INDELs detected',
            'min': 0,
            'scale': 'BrBg',
            'format': '{:,.0f}',
            'shared_key': 'total_variants'
        }
        headers['bases_failed_strand_filter'] = {
            'title': 'Fail Strand',
            'description': 'Total number variants failing the strand-filter.',
            'min': 0,
            'scale': 'YlOrBr',
            'format': '{:,.0f}'
        }
        headers['variant_pos_reported_total'] = {
            'title': 'Total variants reported',
            'description': 'Total number of variants reported.',
            'min': 0,
            'scale': 'Spectral',
            'format': '{:,.0f}',
            'shared_key': 'reported_variants'
        }
        headers['variant_pos_reported_snps'] = {
            'title': 'SNPs reported',
            'description': 'Total number of SNPs reported.',
            'min': 0,
            'scale': 'Spectral',
            'format': '{:,.0f}',
            'shared_key': 'reported_variants'
        }
        headers['variant_pos_reported_indels'] = {
            'title': 'INDELs reported',
            'description': 'Total number INDELs reported.',
            'min': 0,
            'scale': 'Spectral',
            'format': '{:,.0f}',
            'shared_key': 'reported_variants'
        }

        self.general_stats_addcols(self.varscan2_data, headers)
    
    def varscan2_report_plot (self):

        """ Make the HighCharts HTML to plot the reported SNPs"""
        #146 variant positions (106 SNP, 40 indel)
        #12 were failed by the strand-filter
        #99 variant positions reported (99 SNP, 0 indel)

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['variant_pos_snps'] =   { 'name': 'Variants (SNPs)' }
        keys['variant_pos_indels'] = { 'name': 'Variants (INDELs)' }

        # Config for the plot
        config = {
            'id': 'varscan2_total_variants',
            'title': 'VarScan2: Total Variants detected',
            'ylab': '# Variants',
            'cpswitch_counts_label': 'Number of Variants',
            'hide_zero_cats': False
        }

        return bargraph.plot(self.varscan2_data, keys, config)
    
    def varscan2_reported_plot (self):
        """ Make the HighCharts HTML to plot the reported SNPs"""

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['variant_pos_reported_snps'] =   { 'name': 'Variants reported (SNPs)', 'color': '#8bbc21' }
        keys['variant_pos_reported_indels'] = { 'name': 'Variants reported (INDELs)', 'color': '#f7a35c' }

        # Config for the plot
        config = {
            'id': 'varscan2_total_variants_reported',
            'title': 'VarScan2: Total Variants reported',
            'ylab': '# Variants',
            'cpswitch_counts_label': 'Number of Variants',
            'hide_zero_cats': False
        }

        return bargraph.plot(self.varscan2_data, keys, config)
