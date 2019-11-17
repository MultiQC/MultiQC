#!/usr/bin/env python

""" MultiQC module to parse output from MultiVCFAnalyzer """

from __future__ import print_function
from collections import OrderedDict
import logging
import json

from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ MultiVCFAnalyzer module """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='MultiVCFAnalyzer', anchor='multivcfanalyzer',
        href="https://github.com/alexherbig/MultiVCFAnalyzer",
        info="""MultiVCFAnalyzer reads multiple VCF files as produced by the GATK UnifiedGenotyper and after filtering provides the combined genotype calls in a number of formats that are suitable for follow-up analyses such as phylogenetic reconstruction, SNP effect analyses, population genetic analyses etc.""")

        # Find and load any MultiVCFAnalyzer reports
        self.mvcf_data = dict()

        # Find and load JSON file
        for f in self.find_log_files('multivcfanalyzer', filehandles=True):
            self.parse_data(f)

        # Filter samples
        self.mvcf_data = self.ignore_samples(self.mvcf_data)

        # Return if no samples found
        if len(self.mvcf_data) == 0:
            raise UserWarning

        # Add in extra columns to data file
        self.compute_perc_hets()

        # Save data output file
        self.write_data_file(self.mvcf_data, 'multiqc_multivcfanalyzer')

        # Add to General Statistics
        self.addSummaryMetrics()


    def parse_data(self, f):
        try:
            data = json.load(f['f'])
        except Exception as e:
            log.debug(e)
            log.warn("Could not parse MultiVCFAnalyzer JSON: '{}'".format(f['fn']))
            return

        # Parse JSON data to a dict
        for s_name, metrics in data.get('metrics', {}).items():
            if (s_name == 'metadata'):
                continue
            
            if (s_name == 'metrics'):
                for sample in data['metrics'].items():
                    s_clean = sample[0]
                    if s_clean in self.mvcf_data:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_clean))

                    self.add_data_source(f, s_clean)
                    self.mvcf_data[s_clean] = dict()

                    for k, v in sample[1].items():
                        try:
                            self.mvcf_data[s_clean][k] = float(v)
                        except ValueError:
                            self.mvcf_data[s_clean][k] = v

    #Compute % heterozygous snp alleles and add to data
    def compute_perc_hets(self):
        """Take the parsed stats from MultiVCFAnalyzer and add one column per sample """
        for sample in self.mvcf_data:
            try:
                self.mvcf_data[sample]['Heterozygous SNP alleles (percent)'] = (self.mvcf_data[sample]['SNP Calls (het)'] / self.mvcf_data[sample]['SNP Calls (all)'])*100
            except ZeroDivisionError:
                self.mvcf_data[sample]['Heterozygous SNP alleles (percent)'] = 'NA'

    def addSummaryMetrics(self):
        """ Take the parsed stats from MultiVCFAnalyzer and add it to the main plot """

        headers = OrderedDict()
        headers['SNP Calls (all)'] = {
            'title': 'SNPs',
            'description': 'Total number of non-reference calls made',
            'scale': 'OrRd',
            'shared_key': 'snp_call'
        }
        headers['SNP Calls (het)'] = {
            'title': 'Het SNPs',
            'description': 'Total number of non-reference calls not passing homozygosity thresholds',
            'scale': 'OrRd',
            'hidden': True,
            'shared_key': 'snp_call'
        }
        headers['% Hets'] = {
            'title': 'Heterozygous SNP alleles (percent)',
            'description': 'Percentage of heterozygous SNP alleles',
            'scale': 'OrRd',
            'shared_key': 'snp_call'
        }
        headers['allPos'] = {
            'title': 'Bases in SNP Alignment',
            'description': 'Length of FASTA file in base pairs (bp)',
            'scale': 'BuPu',
            'hidden': True,
            'shared_key': 'calls'
        }
        headers['discardedVarCall'] = {
            'title': 'Discarded SNP Call',
            'description': 'Number of non-reference positions not reaching genotyping quality threshold',
            'scale': 'BuPu',
            'hidden': True,
            'shared_key': 'calls'
        }
        headers['filteredVarCall'] = {
            'title': 'Filtered SNP Call',
            'description': 'Number of positions ignored defined in user-supplied filter list',
            'scale': 'BuPu',
            'hidden': True,
            'shared_key': 'calls'
        }
        headers['refCall'] = {
            'title': 'Number of Reference Calls',
            'description': 'Number of reference calls made',
            'scale': 'BuPu',
            'hidden': True,
            'shared_key': 'calls'
        }
        headers['discardedRefCall'] = {
            'title': 'Discarded Reference Call',
            'description': 'Number of reference positions not reaching genotyping quality threshold',
            'scale': 'BuPu',
            'hidden': True,
            'shared_key': 'calls'
        }
        headers['noCall'] = {
            'title': 'Positions with No Call',
            'description': 'Number of positions with no call made as reported by GATK',
            'scale': 'BuPu',
            'hidden': True,
            'shared_key': 'calls'
        }
        headers['coverage (fold)'] = {
            'title': 'SNP Coverage',
            'description': 'Average number of reads covering final calls',
            'scale': 'OrRd',
            'shared_key': 'coverage'
        }
        headers['coverage (percent)'] = {
            'title': '% SNPs Covered',
            'description': 'Percent coverage of all positions with final calls',
            'scale': 'PuBuGn',
            'hidden': True,
            'shared_key': 'coverage'
        }
        headers['unhandledGenotype'] = {
            'title': 'Unhandled Genotypes',
            'description': 'Number of positions discarded due to presence of more than one alternate allele',
            'scale': 'BuPu',
            'hidden': True,
            'shared_key': 'snp_count'
        }

        self.general_stats_addcols(self.mvcf_data, headers)
