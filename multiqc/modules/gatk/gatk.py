#!/usr/bin/env python
""" MultiQC module to parse output from GATK """
from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

# Import the GATK submodules
# import varianteval
from .varianteval import VariantEvalMixin
from .base_recalibrator import BaseRecalibratorMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, BaseRecalibratorMixin, VariantEvalMixin):
    """ GATK has a number of different commands and outputs.
    This MultiQC module supports some but not all. The code for
    each script is split into its own file and adds a section to
    the module output if logs are found. """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='GATK', anchor='gatk', target='GATK',
            href='https://www.broadinstitute.org/gatk/',
            info=(" is a toolkit offering a wide variety of tools with a "
                  "primary focus on variant discovery and genotyping."))

        # Set up class objects to hold parsed data
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()

        # Call submodule functions
        n_reports_found = 0
        n_reports_found += self.parse_gatk_base_recalibrator()
        n_reports_found += self.parse_gatk_varianteval()

        # Exit if we didn't find anything
        if n_reports_found == 0:
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)

    def parse_report(self, lines, table_names):
        """ Parse a GATK report https://software.broadinstitute.org/gatk/documentation/article.php?id=1244

        Only GATTable entries are parsed.  Tables are returned as a dict of tables.
        Each table is a dict of arrays, where names correspond to column names, and arrays
        correspond to column values.

        Args:
            lines (file handle): an iterable over the lines of a GATK report.
            table_names (dict): a dict with keys that are GATK report table names
                (e.g. "#:GATKTable:Quantized:Quality quantization map"), and values that are the
                keys in the returned dict.

        Returns:
            {
                table_1:
                    {
                        col_1: [ val_1, val_2, ... ]
                        col_2: [ val_1, val_2, ... ]
                        ...
                    }
                table_2:
                    ...
            }
        """

        report = dict()

        lines = (l for l in lines)
        for line in lines:
            line = line.rstrip()
            if line in table_names.keys():
                report[table_names[line]] = self.parse_gatk_report_table(lines)
        return report

    def parse_gatk_report_table(self, lines):
        headers = next(lines).rstrip().split()
        table = OrderedDict([(h, []) for h in headers])
        for line in lines:
            line = line.rstrip()

            # testing to see if we have reached the end of a table in a GATKReport
            if line == '':
                break

            for index, value in enumerate(line.split()):
                table[headers[index]].append(value)
        return table
