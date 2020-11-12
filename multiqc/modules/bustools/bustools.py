#!/usr/bin/env python

""" MultiQC module to parse output from bustools inspect """

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import json

from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Bustools', anchor='bustools',
                                            href='https://bustools.github.io/',
                                            info='BUS format is a file format for single-cell RNA-seq data '
                                                 'designed to facilitate the development of modular workflows for data '
                                                 'processing. It consists of a binary representation of barcode and '
                                                 'UMI sequences from scRNA-seq reads, along with sets of equivalence '
                                                 'classes obtained by pseudoalignment of reads to a reference '
                                                 'transcriptome (hence the acronym Barcode, UMI, Set). BUS files are '
                                                 'a convenient and useful checkpoint during single-cell RNA-seq '
                                                 'processing.')

        # Parse logs
        self.bustools_data = dict()
        for f in self.find_log_files('bustools', filehandles=True):
            content = json.load(f['f'])
            s_name = self.clean_s_name(os.path.basename(f['root']), os.path.dirname(f['root']))
            self.bustools_data[s_name] = content

        # Filter to strip out ignored sample names
        self.bustools_data = self.ignore_samples(self.bustools_data)

        if len(self.bustools_data) == 0:
            raise UserWarning

        log.info('Found {} logs'.format(len(self.bustools_data)))
        self.write_data_file(self.bustools_data, 'multiqc_macs')

        self.bustools_general_stats()
        self.bustools_section()

    def bustools_general_stats(self):
        """ Add columns to General Statistics table """
        headers = OrderedDict()
        headers['numRecords'] = {
            'title': 'Bus Records',
            'description': 'The number of Bus Records',
            'min': 0,
            'hidden': False,
            'format': '{:,.0f}'
        }
        headers['numReads'] = {
            'title': 'Reads',
            'description': 'Total number of reads',
            'min': 0,
            'hidden': True,
            'format': '{:,.0f}',
        }
        headers['numBarcodes'] = {
            'title': 'barcodes',
            'description': 'Number of distinct barcodes',
            'min': 0,
            'hidden': False,
            'format': '{:,.0f}',
        }
        headers['meanReadsPerBarcode'] = {
            'title': 'Mean reads per barcode',
            'min': 0,
            'hidden': False,
            'format': '{:,.2f}',
        }
        headers['numUMIs'] = {
            'title': 'distinct UMIs',
            'description': 'Number of distinct Unique Molecular Identifiers (UMIs)',
            'min': 0,
            'hidden': False,
            'format': '{:,.0f}',
        }
        headers['numBarcodeUMIs'] = {
            'title': 'distinct barcode-UMI',
            'description': 'Number of distinct barcode and Unique Molecular Identifiers (UMIs) pairs',
            'min': 0,
            'hidden': False,
            'format': '{:,.0f}',
        }
        headers['meanUMIsPerBarcode'] = {
            'title': 'Mean UMIs per barcode',
            'min': 0,
            'hidden': False,
            'format': '{:,.2f}',
        }
        headers['gtRecords'] = {
            'title': '2xdepth records',
            'description': 'Estimated number of new records at 2x sequencing depth',
            'min': 0,
            'hidden': False,
            'format': '{:,.0f}',
        }
        headers['numTargets'] = {
            'title': 'distinct targets',
            'description': 'Number of distinct targets detected',
            'min': 0,
            'hidden': True,
            'format': '{:,.0f}',
        }
        headers['meanTargetsPerSet'] = {
            'title': 'mean targets',
            'description': 'Mean number of targets per set',
            'min': 0,
            'hidden': True,
            'format': '{:,.2f}',
        }
        headers['numSingleton'] = {
            'title': 'singleton target',
            'description': 'Number of reads with singleton target',
            'min': 0,
            'hidden': True,
            'format': '{:,.2f}',
        }
        headers['gtTargets'] = {
            'title': '2xdepth targets',
            'description': 'Estimated number of new targets at 2x sequencing depth',
            'min': 0,
            'hidden': True,
            'format': '{:,.2f}',
        }
        headers['numBarcodesOnWhitelist'] = {
            'title': 'Whitelisted barcodes',
            'description': 'Number of barcodes in agreement with whitelist',
            'min': 0,
            'hidden': True,
            'format': '{:,.0f}',
        }
        headers['percentageBarcodesOnWhitelist'] = {
            'title': 'Perc. whitelisted barcodes',
            'min': 0,
            'max': 100,
            'scale': 'RdYlGn',
            'hidden': False,
            'format': '{:,.2f}',
        }
        headers['numReadsOnWhitelist'] = {
            'title': 'Whitelisted reads',
            'description': 'Number of reads with barcode in agreement with whitelist',
            'min': 0,
            'hidden': True,
            'format': '{:,.0f}',
        }
        headers['percentageReadsOnWhitelist'] = {
            'title': 'Perc. whitelisted reads',
            'min': 0,
            'max': 100,
            'scale': 'RdYlGn',
            'hidden': True,
            'format': '{:,.2f}',
        }
        self.general_stats_addcols(self.bustools_data, headers)

    def bustools_section(self):
        """ Add bargraphs showing the mean UMIs per barcode and percentages in whitelist """
        mean_umis = {sample: {'UMIs per barcode': values['meanUMIsPerBarcode']} for sample, values in
                     self.bustools_data.items()}

        percentage_whitelist = {sample: {'Reads on whitelist': values['percentageReadsOnWhitelist'],
                                         'Barcodes on whitelist': values['percentageBarcodesOnWhitelist']}
                                for sample, values in self.bustools_data.items()}

        self.add_section(
            name='Mean number of UMIs per barcode',
            anchor='bustools-umis',
            description='Each unique barcode represents a cell and each Unique Molecular Identifier (UMI) represents '
                        'a unique transcript molecule. By counting the mean number of UMIs per barcode, you '
                        'effectively calculate the average number of unique transcripts per cell.',
            plot=bargraph.plot(mean_umis,
                               pconfig={
                                        'id': 'bus_umis',
                                        'title': 'Bustools: Mean number of UMIs per barcode per sample',
                                        'cpswitch': False,
                                        'tt_percentages': False,
                                        'ylab': 'Mean UMIs per barcode'
                                        }
                               )
        )

        self.add_section(
            name='Percentage in whitelist',
            anchor='bustools-reads',
            description='The whitelist is a list of unique barcodes used in your protocol, either provided or inferred '
                        'from the data. Each unique barcode from the whitelist represents a cell. The percentage of '
                        'reads with barcode / barcodes in the whitelist is a measure of percentage of reads that could '
                        'be asigned to a cell.',
            plot=bargraph.plot(percentage_whitelist,
                               pconfig={'id': 'bus_reads',
                                        'title': 'Bustools: percentage of barcodes and/or reads with barcodes in the '
                                                 'whitelist',
                                        'ymax': 100,
                                        'ymix': 0,
                                        'cpswitch': False,
                                        'tt_percentages': False,
                                        'ylab': 'Percentage of reads with barcodes in the whitelist',
                                        'stacking': None,
                                        'data_labels': ['Reads', 'Bases']
                                        }
                               )
        )
