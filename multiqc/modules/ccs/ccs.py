#!/usr/bin/env python

""" MultiQC module to parse output from CCS """

import logging
import re

from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialse the parent object
        super(MultiqcModule, self).__init__(name='ccs', anchor='ccs',
                href='https://github.com/PacificBiosciences/ccs',
                info='Generate highly accurate single-molecule consensus '
                     'reads')

        for f in self.find_log_files('ccs', filehandles=True):
            parse_ccs_log(f['f'])


def parse_ccs_log(file_content):
    """ Parse ccs log file """
    data = dict()
    key, value = parse_line(next(file_content))
    data[key] = value
    return data


def parse_line(line):
    """ Parse a line from the ccs log file """
    # Split the line on the colon character
    key, value = line.strip().split(':')

    # Remove (A), (B), etc annotations from the key field
    key = re.sub(r' + [(]A[)] +','', key)

    # Remove the percentage between bracets from the value
    value = re.sub(r' + [(]\d+.*[)]','', value.strip())

    # All values are integers
    value = int(value)

    return key, value
