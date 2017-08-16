#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" MultiQC submodule to parse output from GATK BaseRecalibrator """

import logging
from collections import OrderedDict, defaultdict
from multiqc.plots import bargraph, table, linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class BaseRecalibratorMixin():
    def parse_gatk_base_recalibrator(self):
        """ Find GATK BaseRecalibrator logs and parse their data """

        self.gatk_base_recalibrator = dict()
        for f in self.find_log_files('gatk/base_recalibrator', filehandles=True):
            parsed_data = self.parse_report(
                f['f'],
                {
                    '#:GATKTable:Arguments:Recalibration argument collection values used in this run': 'arguments',
                    '#:GATKTable:Quantized:Quality quantization map': 'quality_quantization_map',
                    '#:GATKTable:RecalTable0:': 'recal_table_0',
                    '#:GATKTable:RecalTable1:': 'recal_table_1',
                    '#:GATKTable:RecalTable2:': 'recal_table_2',
                }
            )
            if len(parsed_data) > 0:
                if f['s_name'] in self.gatk_base_recalibrator:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
                self.add_data_source(f, section='base_recalibrator')
                self.gatk_base_recalibrator[f['s_name']] = parsed_data

        # Filter to strip out ignored sample names
        self.gatk_base_recalibrator = self.ignore_samples(self.gatk_base_recalibrator)

        if len(self.gatk_base_recalibrator) > 0:
            # Write parsed report data to a file (restructure first)
            self.write_data_file(self.gatk_base_recalibrator, 'multiqc_gatk_base_recalibrator')

            # Reported vs empirical quality scores
            self.add_section(
                name='Reported vs Empirical Quality Scores',
                anchor='gatk-reported-vs-empirical-quality-scores',
                plot=quality_score_vs_no_of_observations(self.gatk_base_recalibrator)
            )

        # Return the number of logs that were found
        return len(self.gatk_base_recalibrator)


def quality_score_vs_no_of_observations(data):
    """ Return HTML for the quality score vs number of observations line plot """

    table = data['quality_quantization_map']
    xy = {x: y for x, y in zip(table['QualityScore'], table['Count'])}
    return linegraph.plot({'sample_1': xy})
