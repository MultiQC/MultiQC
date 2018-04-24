#!/usr/bin/env python

""" MultiQC module to parse output from KAT """
import logging
from collections import OrderedDict
import json

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='K-mer Analysis Toolkit', anchor='kat',
                                            href="https://github.com/TGAC/KAT",
                                            info="is an toolkit for analysing sequencing data via its k-mer spectra.")

        # Find and load any KAT dist analysis reports
        self.kat_data = dict()
        for c_file in self.find_log_files('kat'):
            s_name = self.clean_s_name(c_file['s_name'].replace(".dist_analysis", ""), c_file['root'])
            content = json.loads(c_file['f'])
            self.kat_data[s_name] = self.parse_kat_report(content)

        # Filter to strip out ignored sample names
        self.kat_data = self.ignore_samples(self.kat_data)

        if len(self.kat_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.kat_data)))

        # Write parsed report data to a file
        self.write_data_file(self.kat_data, 'multiqc_kat')

        headers = OrderedDict()
        headers["kmer_peaks"] = {
            'title': '# of Kmer Peaks',
            'description': "Number of peaks identified in the K-mer spectra",
            'scale': False,
            'format': '{:,.0f}'
        }
        headers["gc_peaks"] = {
            'title': '# of GC Peaks',
            'description': "Number of peaks identified in the GC distribution",
            'scale': False,
            'format': '{:,.0f}'
        }
        headers["est_genome_size"] = {
            'title': 'Est. genome Size',
            'description': "Estimated Genome Size based on K-mer spectra",
            'scale': "BuPu",
            'format': '{:,.0f}'
        }
        headers["mean_kmer_freq"] = {
            'title': 'Mean K-mer Freq.',
            'description': "Mean K-mer Frequency, provides an estimate of sequencing coverage",
            'scale': "Greens",
            'format': '{:,.0f}',
            'suffix': 'x'
        }

        kat_config = {
            'namespace': 'KAT',
        }

        # Basic Stats Table
        self.add_section(
            name='KAT Distribution Analysis',
            anchor='kat-first',
            description='Table showing k-mer coverage distributions and if available GC distributions',
            helptext="This table can give a quick idea of potential contaminants that can be identified via unexpected numbers of k-mer or gc peaks in the data",
            plot=table.plot(self.kat_data, headers, kat_config)
        )

    def parse_kat_report(self, content):
        table_data = {}
        if "gc" in content and "coverage" in content:
            # GCP
            table_data['kmer_peaks'] = content['coverage']['nb_peaks']
            table_data['mean_kmer_freq'] = content['coverage']['mean_freq']
            table_data['est_genome_size'] = content['coverage']['est_genome_size']
            table_data['gc_peaks'] = content['gc']['nb_peaks']
        elif "main_dist" in content:
            # Spectra CN
            table_data['kmer_peaks'] = content['main_dist']['nb_peaks']
            table_data['mean_kmer_freq'] = content['main_dist']['mean_freq']
            table_data['est_genome_size'] = content['main_dist']['est_genome_size']
            table_data['gc_peaks'] = 0
        elif "k" in content:
            # Hist
            table_data['kmer_peaks'] = content['nb_peaks']
            table_data['mean_kmer_freq'] = content['mean_freq']
            table_data['est_genome_size'] = content['est_genome_size']
            table_data['gc_peaks'] = 0
        else:
            log.error("Unexpected JSON configuration")

        return table_data
