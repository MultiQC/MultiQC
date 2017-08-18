#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MultiQC submodule to parse output from GATK BaseRecalibrator """

import logging
from collections import OrderedDict, defaultdict
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class BaseRecalibratorMixin():
    def parse_gatk_base_recalibrator(self):
        """ Find GATK BaseRecalibrator logs and parse their data """

        report_table_headers = {
            '#:GATKTable:Arguments:Recalibration argument collection values used in this run': 'arguments',
            '#:GATKTable:Quantized:Quality quantization map': 'quality_quantization_map',
            '#:GATKTable:RecalTable0:': 'recal_table_0',
            '#:GATKTable:RecalTable1:': 'recal_table_1',
            '#:GATKTable:RecalTable2:': 'recal_table_2',
        }
        samples_kept = set()
        self.gatk_base_recalibrator = {table_name: {} for table_name in
                                       report_table_headers.values()}
        for f in self.find_log_files('gatk/base_recalibrator', filehandles=True):
            parsed_data = self.parse_report(f['f'].readlines(), report_table_headers)
            if len(parsed_data) > 0:
                if f['s_name'] in samples_kept:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                samples_kept.add(f['s_name'])

                self.add_data_source(f, section='base_recalibrator')
                for table_name, sample_tables in parsed_data.items():
                    self.gatk_base_recalibrator[table_name][f['s_name']] = sample_tables

        # Filter to strip out ignored sample names
        self.gatk_base_recalibrator = {
            table_name: self.ignore_samples(sample_tables)
            for table_name, sample_tables
            in self.gatk_base_recalibrator.items()
        }

        n_reports_found = len(samples_kept)
        if n_reports_found > 0:
            log.info("Found {} BaseRecalibrator reports".format(n_reports_found))

            self.add_quality_score_vs_no_of_observations_section()

        return n_reports_found

    def add_quality_score_vs_no_of_observations_section(self):
        """ Add a section for the quality score vs number of observations line plot """

        sample_tables = self.gatk_base_recalibrator['quality_quantization_map']
        sample_data = {
            sample: {int(x): int(y) for x, y in zip(table['QualityScore'], table['Count'])}
            for sample, table in sample_tables.items()
        }

        plot = linegraph.plot(
            sample_data,
            pconfig={
                'title': "Observed Quality Score Counts",
                'id': 'gatk-base-recalibrator-quality-score-vs-number-of-observations',
                'xlab': 'Observed Quality Score',
                'ylab': 'Count',
                'yDecimals': False,
                'xDecimals': False,
                'tt_label': '{point.x}: {point.y:.0f}'
            })

        # Reported vs empirical quality scores
        self.add_section(
            name='Observed Quality Scores',
            description=(
                'This plot shows the distribution of base quality scores in each sample before and '
                'after base quality score recalibration (BQSR). Applying BQSR should broaden the '
                'distribution of base quality scores.'
            ),
            helptext=(
                'For more information see <a href=https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr>'
                'the Broad\'s description of BQSR</a>.'
            ),
            plot=plot,
        )
