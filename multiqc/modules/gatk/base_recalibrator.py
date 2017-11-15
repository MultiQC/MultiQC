#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MultiQC submodule to parse output from GATK BaseRecalibrator """

import logging
from enum import IntEnum
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class RecalTableType(IntEnum):
    pre_recalibration = 0
    post_recalibration = 1

    def human_readable(self):
        return self.name.capitalize().replace('_', '-')


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
        samples_kept = {rt_type: set() for rt_type in RecalTableType}
        self.gatk_base_recalibrator = {recal_type:
                                           {table_name: {}
                                            for table_name
                                            in report_table_headers.values()}
                                       for recal_type in RecalTableType}

        for f in self.find_log_files('gatk/base_recalibrator', filehandles=True):
            parsed_data = self.parse_report(f['f'].readlines(), report_table_headers)
            rt_type = determine_recal_table_type(parsed_data)
            if len(parsed_data) > 0:
                if f['s_name'] in samples_kept[rt_type]:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                else:
                    samples_kept[rt_type].add(f['s_name'])

                self.add_data_source(f, section='base_recalibrator')
                for table_name, sample_tables in parsed_data.items():
                    self.gatk_base_recalibrator[rt_type][table_name][
                        f['s_name']] = sample_tables

        # Filter to strip out ignored sample names
        for rt_type in RecalTableType:
            for table_name, sample_tables in self.gatk_base_recalibrator[rt_type].items():
                self.gatk_base_recalibrator[rt_type][table_name] = self.ignore_samples(
                    sample_tables)

        n_reports_found = sum([len(samples_kept[rt_type]) for rt_type in RecalTableType])

        if n_reports_found > 0:
            log.info("Found {} BaseRecalibrator reports".format(n_reports_found))

            self.add_quality_score_vs_no_of_observations_section()

        return n_reports_found

    def add_quality_score_vs_no_of_observations_section(self):
        """ Add a section for the quality score vs number of observations line plot """

        sample_data = []
        data_labels = []
        for rt_type in RecalTableType:
            sample_tables = self.gatk_base_recalibrator[rt_type]['quality_quantization_map']
            if len(sample_tables) == 0:
                continue

            sample_data.append({
                sample: {int(x): int(y) for x, y in zip(table['QualityScore'], table['Count'])}
                for sample, table in sample_tables.items()
            })

            sample_y_sums = {
                sample: sum(int(y) for y in table['Count'])
                for sample, table
                in sample_tables.items()
            }

            sample_data.append({
                sample: {
                    int(x): float(y) / sample_y_sums[sample]
                    for x, y in zip(table['QualityScore'], table['Count'])
                }
                for sample, table in sample_tables.items()
            })

            flat_proportions = [float(y) / sample_y_sums[sample]
                                for sample, table in sample_tables.items()
                                for y in table['Count']]
            prop_ymax = max(flat_proportions)
            data_labels.append({'name': "{} Count".format(rt_type.human_readable()),
                                'ylab': 'Count'})
            data_labels.append({'ymax': prop_ymax,
                                'name': "{} Percent".format(rt_type.human_readable()),
                                'ylab': 'Percent'})

        plot = linegraph.plot(
            sample_data,
            pconfig={
                'title': "Observed Quality Score Counts",
                'id': 'gatk-base-recalibrator-quality-score-vs-number-of-observations',
                'xlab': 'Observed Quality Score',
                'ylab': 'Count',
                'xDecimals': False,
                'data_labels': data_labels,
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
                'For more information see '
                '[the Broad\'s description of BQSR]'
                '(https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr)'
                '.'
            ),
            plot=plot,
        )


def determine_recal_table_type(parsed_data):
    arguments = dict(zip(parsed_data['arguments']['Argument'], parsed_data['arguments']['Value']))
    if arguments['recalibration_report'] == 'null':
        return RecalTableType.pre_recalibration
    else:
        return RecalTableType.post_recalibration
