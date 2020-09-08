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


def parse_PacBio_log(file_content):
    """ Parse ccs log file """
    data = dict()
    key, value = parse_line(next(file_content))
    data[key] = value
    return data


def parse_line(line):
    """ Parse a line from the ccs log file """
    data = dict()

    # If we got an empty line to parse
    if not line.strip():
        return data

    # Split the line on the colon character
    key, values = line.strip().split(':')

    # The key can have multiple parts
    keys = key.strip().split()

    # Does the key have an annotation (A), (B) etc
    if re.fullmatch('[(][A-Z][)]', keys[-1]):
        # We store the annotation without the bracets
        data['annotation'] = keys[-1][1:-1]
        # And we add the rest of the key as name
        data['name'] = ' '.join(keys[:-1])
    # Otherwise, we just store the name
    else:
        data['name'] = ' '.join(keys)

    return data
