#!/usr/bin/env python
""" MultiQC module to parse output from GATK """
from __future__ import print_function
from collections import OrderedDict, defaultdict
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
        n = dict()

        # Call submodule functions
        n['varianteval'] = self.parse_gatk_varianteval()
        if n['varianteval'] > 0:
            log.info("Found {} VariantEval reports".format(n['varianteval']))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)

    def parse_report(self, f, table_names):
        """ Parse a GATK report https://software.broadinstitute.org/gatk/documentation/article.php?id=1244

        Only GATTable entries are parsed.  Tables are returned as a dict of tables.
        Each table is a dict of arrays, where names correspond to column names, and arrays
        correspond to column values.

        Args:
            f (file handle): a file handle to a GATK report.
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

        data = dict()
        while True:
            line = f.readline()
            if line == '':
                break
            line = line.rstrip()
            if line in table_names.keys():
                data[table_names[line]] = self.parse_gatk_report_table(f)
        return data

    def parse_gatk_report_table(self, f):
        data = defaultdict(list)
        headers = f.readline().rstrip().split()
        while True:
            line = f.readline()
            line = line.rstrip()
            if line == '':
                break
            for index, value in enumerate(line.split()):
                data[headers[index]].append(value)
        return data
