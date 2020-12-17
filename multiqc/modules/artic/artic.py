#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MultiQC module to parse output files from ARTIC """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

# Alignment_Length_Threshold drops binned reads that are <X% of amplicon length)
Alignment_Length_Threshold = 0.95

# Amplicon_Dropout_Val will report amplicon dropout in any amplicon which has fewer than X reads
Amplicon_Dropout_Val = 50


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='ARTIC',
            anchor='ARTIC',
            href='https://github.com/artic-network/fieldbioinformatics',
            info="is a bioinformatics pipeline for working with virus sequencing data sequenced with nanopore."
        )

        # Find and load ARTIC align_trim results
        self.aligntrim_data = dict()
        self.vcfcheck_data = dict()
        self.artic_stats = dict()
        for f in self.find_log_files('artic/aligntrim_reports', filehandles=True):
            self.parse_aligntrim_report(f)
        for f in self.find_log_files('artic/vcf_reports', filehandles=True):
            self.parse_vcf_report(f)

        # Filter to strip out ignored sample names
        self.aligntrim_data = self.ignore_samples(self.aligntrim_data)
        self.vcfcheck_data = self.ignore_samples(self.vcfcheck_data)

        # Stop if no files are found
        num_samples = max(len(self.aligntrim_data), len(self.vcfcheck_data))
        if num_samples == 0:
            raise UserWarning

        # Log number of reports
        log.info("Found {} align_trim reports".format(len(self.aligntrim_data)))
        log.info("Found {} vcf_check reports".format(len(self.vcfcheck_data)))

        # Write parsed data to a file
        self.write_data_file(self.aligntrim_data,
                             'multiqc_artic_aligntrim_summary')
        self.write_data_file(self.vcfcheck_data,
                             'multiqc_artic_vcfcheck_summary')

        # Process collected files
        self.process_artic_data()

        # Update General Stats Table
        self.artic_general_stats_table()

        # Add amplicon coverage lot
        self.plot_amplicons()

    # Parse an aligntrim report
    def parse_aligntrim_report(self, f):
        """ Parse ARTIC aligntrim report file
        """
        sample = f['s_name'].replace('.alignreport', '')

        # Check if overwriting
        if sample in self.aligntrim_data:
            log.debug(
                "Duplicate sample name found! Overwriting: {}".format(sample))

        # Set up data holder for sample
        self.aligntrim_data[sample] = dict()

        # Skip first line of report (header)
        _ = f['f'].readline()

        # Read aligntrim report and bin reads by amplicon
        for l in f['f']:
            fields = l.rstrip().split('\t')

            # Check read is from a properly paired amplicon
            if int(fields[12]) != 1:
                continue

            # Check the read alignment length covers enough of the amplicon
            aLen = int(fields[11]) - int(fields[10])
            rLen = int(fields[2]) - int(fields[1])
            if aLen < (Alignment_Length_Threshold * rLen):
                continue

            # Update the read count for this amplicon
            if fields[3] not in self.aligntrim_data[sample]:
                self.aligntrim_data[sample][fields[3]] = 0
            self.aligntrim_data[sample][fields[3]] += 1

    # Parse a vcf_check report
    def parse_vcf_report(self, f):
        """ Parse ARTIC vcf_check report file
        """
        sample = f['s_name']

        # Check if overwriting
        if sample in self.vcfcheck_data:
            log.debug("Duplicate sample name found! Overwriting: {}".format(sample))

        # Set up data holder for sample
        self.vcfcheck_data[sample] = dict()

        # Read vcfcheck report and get important stuff out (NOTE: more to be added in next release)
        total_vars = 0
        passed_vars = 0
        for l in f['f']:
            match = re.search(r'.*\t(\d+)\svariant\srecords\sprocessed', l)
            if match:
                total_vars = int(match.group(1))
            match = re.search(r'.*\t(\d+)\svariant\srecords\spassed\schecks', l)
            if match:
                passed_vars = int(match.group(1))
        self.vcfcheck_data[sample]['overlap_fails'] = total_vars - passed_vars

    # Process collected data for all samples
    def process_artic_data(self):
        """ Process the data from the parsed artic files
        """

        # Get the amplicon counts
        for sample in self.aligntrim_data:
            self.artic_stats[sample] = dict()
            self.artic_stats[sample]['amplicon_dropouts'] = 0
            for count in self.aligntrim_data[sample].values():
                if count < Amplicon_Dropout_Val:
                    self.artic_stats[sample]['amplicon_dropouts'] += 1

        # Get the vcf counts
        for sample in self.vcfcheck_data:

            # check for that sample already has an aligntrim report
            if sample not in self.artic_stats:
                log.debug("Sample has vcf report but no aligntrim report! {}".format(sample))
                self.artic_stats[sample] = dict()
            self.artic_stats[sample]['overlap_fails'] = self.vcfcheck_data[sample]['overlap_fails']

    # Add to general stats table
    def artic_general_stats_table(self):
        """ Take the parsed stats from the artic reports and add to the
        general stats table"""

        headers = OrderedDict()
        headers['amplicon_dropouts'] = {
            'title': '# low cov. amplicons',
            'description': 'The number of amplicons for this sample with read coverage <{}' .format(Amplicon_Dropout_Val),
            'min': 0,
            'scale': 'RdYlGn-rev',
            'format': '{:,.0f}'
        }
        headers['overlap_fails'] = {
            'title': '# overlap var. fails',
            'description': 'The number of variants detected only once in scheme overlap regions',
            'min': 0,
            'scale': 'RdYlGn-rev',
            'format': '{:,.0f}'
        }
        self.general_stats_addcols(self.artic_stats, headers)

    # Create amplicon coverage plot
    def plot_amplicons(self):
        """ Plot the amplicon coverage"""

        pconfig = {
            'id': 'amplicon_plot',
            'title': 'Module: Amplicon Coverage',
            'ylab': '# reads',
            'xlab': 'amplicon',
            'yPlotLines': [{
                'color': '#FF0000',
                'width': 2,
                'dashStyle': 'LongDash',
                'label': 'Amplicon dropout',
                'value': Amplicon_Dropout_Val
            }],
            'categories': True
        }

        if self.aligntrim_data is not None:
            amplicon_plot_html = linegraph.plot(
                self.aligntrim_data, pconfig)
            self.add_section(
                name='Amplicon coverage',
                anchor='artic_amplicon_plot',
                # description='Mean depth of coverage per amplicon.',
                helptext='''
                This plot summarises the number of reads that were assigned to each amplicon in the primer scheme.

                We use the `align_trim` report file from the ARTIC pipeline and group each read by its assigned amplicon.

                If the length of alignment between read and reference is <{}% of the amplicon length, the read discarded from the coverage plot.

                If the total number of reads assigned to an amplicon is below {} (red dashed line), the amplicon is marked as dropped out.
                ''' .format(Alignment_Length_Threshold*100, Amplicon_Dropout_Val),
                plot=amplicon_plot_html
            )
