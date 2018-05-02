#!/usr/bin/env python

""" MultiQC module to parse output from Longranger """

from __future__ import print_function
from collections import OrderedDict
import logging
import re
import os
import csv

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ Longranger module """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Longranger', anchor='longranger',
        href="https://www.10xgenomics.com/",
        info="A set of analysis pipelines that perform sample demultiplexing, "
        "barcode processing, alignment, quality control, variant calling, phasing, "
        "and structural variant calling.")

        self.headers = OrderedDict()
        self.headers['longranger_version'] = {
                'scale': False,
                'title': 'Version',
                'description': 'The version of the Long Ranger software used to generate the results.',
                'hidden': True
        }
        self.headers['instrument_ids'] = {
                'scale': False,
                'title': 'Instrument ID',
                'description': 'The list of instrument IDs used to generate the input reads.',
                'hidden': True
        }
        self.headers['gems_detected'] = {
                'scale': 'RdYlGn',
                'title': '# gems',
                'modify': lambda x: float(x) / 1000000.0,
                'suffix': 'M',
                'description': 'The number of Chromium GEMs that were collected and which generated a non-trivial number of read-pairs.'
        }
        self.headers['mean_dna_per_gem'] = {
                'scale': 'RdYlGn',
                'title': 'DNA per gem',
                'modify': lambda x: float(x) / 1000000.0,
                'suffix': 'M',
                'description': 'The average number of base pairs of genomic DNA loaded into each GEM. This metric is based on the observed extents of read-pairs on each molecule.',
                'hidden': True
        }
        self.headers['bc_on_whitelist'] = {
                'scale': 'RdYlGn',
                'title': 'Valid BCs',
                'modify': lambda x: float(x) * 100.0,
                'suffix': '%',
                'description': 'The fraction of reads that carried a valid 10x barcode sequence.',
                'hidden': True
        }
        self.headers['bc_mean_qscore'] = {
                'scale': 'RdYlGn',
                'title': 'BC Qscore',
                'description': 'The mean base quality value on the barcode bases.',
                'hidden': True
        }
        self.headers['n50_linked_reads_per_molecule'] = {
                'scale': 'RdYlGn',
                'title': 'N50 read per mol.',
                'description': 'The N50 number of read-pairs per input DNA molecule. Equivalently, half of read-pairs came from molecules with this many or greater read-pairs, and half came from molecules with fewer read pairs.',
                'hidden': True
        }
        self.headers['corrected_loaded_mass_ng'] = {
                'scale': 'RdYlGn',
                'title': 'Loaded',
                'description': 'The estimated number of nanograms of DNA loaded into the input well of the Chromium chip. This metric is calculated by measuring the mean amount of DNA covered by input molecules in each GEM, then multiplying by the ratio of the chip input to the sample volume in each GEM.',
                'suffix': 'ng'
        }
        self.headers['loaded_mass_ng'] = {
                'scale': 'RdYlGn',
                'title': 'Loaded',
                'description': 'This metric was found to overestimate the true loading by a factor of 1.6, due primarily to denaturation of the input DNA.',
                'suffix': 'ng'
        }
        self.headers['genes_phased_lt_100kb'] = {
                'scale': 'RdYlGn',
                'title': '# genes phased < 100kb',
                'description': 'Fraction of genes shorter than 100kb with >1 heterozygous SNP that are phased into a single phase block.',
                'hidden': True,
                'modify': lambda x: float(x) * 100.0,
                'suffix': '%'
        }
        self.headers['longest_phase_block'] = {
                'scale': 'RdYlGn',
                'title': 'Longest phased',
                'modify': lambda x: float(x) / 1000000.0,
                'suffix': 'Mbp',
                'description': 'Size of the longest phase block, in base pairs',
                'hidden': True
        }
        self.headers['n50_phase_block'] = {
                'scale': 'RdYlGn',
                'title': 'N50 phased',
                'modify': lambda x: float(x) / 1000000.0,
                'suffix': 'Mbp',
                'description': 'N50 length of the called phase blocks, in base pairs.'
        }
        self.headers['molecule_length_mean'] = {
                'scale': 'RdYlGn',
                'title': 'Mol size',
                'modify': lambda x: float(x) / 1000.0,
                'suffix': 'Kbp',
                'description': 'The length-weighted mean input DNA length in base pairs.'
        }
        self.headers['molecule_length_stddev'] = {
                'scale': 'RdYlGn',
                'title': 'Mol stddev',
                'modify': lambda x: float(x) / 1000.0,
                'suffix': 'Kbp',
                'description': 'The length-weighted standard deviation of the input DNA length distribution in base pairs.',
                'hidden': True
        }
        self.headers['number_reads'] = {
                'scale': 'PuBu',
                'title': '# Reads',
                'modify': lambda x: float(x) / 1000000.0,
                'suffix': 'M',
                'description': 'Total number of reads supplied to Long Ranger.'
        }
        self.headers['median_insert_size'] = {
                'scale': 'PuBu',
                'title': 'Insert size',
                'suffix': 'bp',
                'description': 'Median insert size of aligned read pairs.',
                'format': '{:,.0f}',
                'hidden': True
        }
        self.headers['mean_depth'] = {
                'scale': 'PuBu',
                'title': 'Depth',
                'description': 'Mean read depth, including PCR duplicate reads. In WGS mode, this is measured across the genome; in targeted mode, this is the measure inside targeted regions.',
                'suffix': 'X'
        }
        self.headers['zero_coverage'] = {
                'scale': 'PuBu-rev',
                'title': 'Zero cov',
                'modify': lambda x: float(x) * 100.0,
                'suffix': '%',
                'description': 'Fraction of non-N bases in the genome with zero coverage.'
        }
        self.headers['mapped_reads'] = {
                'scale': 'PuBu',
                'title': 'Mapped',
                'modify': lambda x: float(x) * 100.0,
                'suffix': '%',
                'description': 'Number of input reads that were mapped.'
        }
        self.headers['pcr_duplication'] = {
                'scale': 'PuBu-rev',
                'title': 'PCR Dup',
                'description': 'Fraction of reads marked as PCR duplicates. To be marked as PCR duplicates, reads must have the same mapping extents on the genome and the same 10x barcode.',
                'suffix': '%',
                'modify': lambda x: float(x) * 100.0
        }
        self.headers['on_target_bases'] = {
                'scale': 'PuBu',
                'title': 'On target',
                'description': 'Fraction of aligned bases mapped with the target regions in targeted mode. Only bases inside the intervals of target BED file are counted.',
                'suffix': '%',
                'modify': lambda x: 0 if x=="" else float(x) * 100.0
        }
        self.headers['r1_q20_bases_fract'] = {
                'scale': 'PuBu',
                'title': 'R1 Q20',
                'description': 'Fraction of bases in R1 with base quality >= 20.',
                'hidden': True,
                'suffix': '%',
                'modify': lambda x: float(x) * 100.0
        }
        self.headers['r1_q30_bases_fracts'] = {
                'scale': 'PuBu',
                'title': 'R1 Q30',
                'description': 'Fraction of bases in R1 with base quality >= 30.',
                'hidden': True,
                'suffix': '%',
                'modify': lambda x: float(x) * 100.0
        }
        self.headers['r2_q20_bases_fract'] = {
                'scale': 'PuBu',
                'title': 'R2 Q20',
                'description': 'Fraction of bases in R2 with base quality >= 20.',
                'hidden': True,
                'suffix': '%',
                'modify': lambda x: float(x) * 100.0
        }
        self.headers['r2_q30_bases_fract'] = {
                'scale': 'PuBu',
                'title': 'R2 Q30',
                'description': 'Fraction of bases in R2 with base quality >= 30.',
                'hidden': True,
                'suffix': '%',
                'modify': lambda x: float(x) * 100.0
        }
        self.headers['si_q20_bases_fract'] = {
                'scale': 'PuBu',
                'title': 'SI Q20',
                'description': 'Fraction of bases in the sample index with base quality >= 20.',
                'hidden': True,
                'suffix': '%',
                'modify': lambda x: float(x) * 100.0
        }
        self.headers['si_q30_bases_fract'] = {
                'scale': 'PuBu',
                'title': 'SI Q30',
                'description': 'Fraction of bases in the sample index with base quality >= 30.',
                'hidden': True,
                'suffix': '%',
                'modify': lambda x: float(x) * 100.0
        }
        self.headers['bc_q20_bases_fract'] = {
                'scale': 'PuBu',
                'title': 'BC Q20',
                'description': 'Fraction of bases in the barcode with base quality >= 20.',
                'hidden': True,
                'suffix': '%',
                'modify': lambda x: float(x) * 100.0
        }
        self.headers['bc_q30_bases_fract'] = {
                'scale': 'PuBu',
                'title': 'BC Q30',
                'description': 'Fraction of bases in the barcode with base quality >= 30.',
                'hidden': True,
                'suffix': '%',
                'modify': lambda x: float(x) * 100.0
        }
        self.headers['snps_phased'] = {
                'scale': 'YlGn',
                'title': 'SNPs phased',
                'description': 'Fraction of called SNPs that were phased.',
                'modify': lambda x: float(x) * 100.0,
                'suffix': '%',
        }
        self.headers['large_sv_calls'] = {
                'scale': 'YlGn',
                'title': 'Large SVs',
                'description': 'Large structural variants called by Longranger. Not including blacklisted regions.',
                'format': '{:,.0f}'
        }
        self.headers['short_deletion_calls'] = {
                'scale': 'YlGn',
                'title': 'Short dels',
                'description': 'Short deletions called by Longranger.',
                'format': '{:,.0f}',
                'hidden': True
        }

        ### Parse the data
        self.longranger_data = dict()
        self.paths_dict = dict()
        for f in self.find_log_files('longranger/invocation'):
            sid = self.parse_invocation(f['f'])
            self.paths_dict[os.path.basename(f['root'])] = sid
        
        running_name = 1
        for f in self.find_log_files('longranger/summary'):
            data = self.parse_summary(f['f'])
            updir, _ = os.path.split(f['root'])
            base_updir = os.path.basename(updir)
            sid = 'longranger#{}'.format(running_name)
            if base_updir in self.paths_dict.keys():
                sid = self.paths_dict[base_updir]
            else:
                log.debug('Did not find _invocation file')
                running_name += 1
            
            self.longranger_data[sid] = data 

        # Filter to strip out ignored sample names
        self.longranger_data = self.ignore_samples(self.longranger_data)

        if len(self.longranger_data) == 0:
            raise UserWarning


        ### Write the report
        self.write_data_file(self.longranger_data, 'multiqc_longranger')
        config_table = {
            'id': 'longranger_table',
            'namespace': 'longranger'
        }
        self.add_section (
            name = 'Longranger statistics',
            anchor = 'longranger-table',
            description = 'Statistics gathered from Longranger reports. ' \
                    'There are more columns available but they are hidden by default.',
            helptext = 'Parses the files `summary.csv` and `_invocation` found in the ' \
                    'output directory of Longranger. If `_invocation` is not found ' \
                    'the sample IDs will be missing and they will be given a running ' \
                    'number. E.g., `longranger#1` and `longranger#2`.',
            plot = table.plot(self.longranger_data, self.headers, config_table)
        )

        log.info("Found {} reports".format(len(self.longranger_data.keys())))

        # Write parsed report data to a file
        self.write_data_file(self.longranger_data, 'multiqc_longranger')


    def parse_invocation(self, content):
        sid_pat = re.compile('    sample_id = \"(.*)\",')
        
        sid = None
        for l in content.splitlines():
            sid_m = re.match(sid_pat,l)
            if sid_m is not None:
                sid = sid_m.groups()[0]
        return sid

    def parse_summary(self, content):
        
        out_dict = OrderedDict()
        lines = content.splitlines()
        data = list(zip(lines[0].strip().split(','), lines[1].strip().split(',')))
        for i,j in data:
            out_dict[i] = j
            
        return out_dict
