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

        self.mod_data = dict()
        for f in self.find_log_files('ccs', filehandles=True):
            self.add_data_source(f)
            self.mod_data[f['s_name']] = parse_PacBio_log(f['f'])
            self.write_data_file(self.mod_data, 'multiqc_ccs')

        # If we found no data
        if not self.mod_data:
            raise UserWarning

        print(self.mod_data)

def parse_PacBio_log(file_content):
    """ Parse ccs log file """
    data = dict()
    # This is a local dictionary to store which annotations belong to which
    # result dictionary. This will be used to structure the output, but will
    # not be part of the output itself
    annotations = dict()
    current_annotation = None

    for line in file_content:
        # Get rid of trailing newlines
        line = line.strip('\n')
        # Did we enter a new section with annotations for an earlier result?
        # If so, we will only add an empty dictionary we the correct name
        section_header_pattern = ' for [(][A-Z][)][:]'
        if re.search(section_header_pattern, line):
            linedata = parse_line(line)
            ann = linedata['annotation']
            # Cut off the ' for (B):' part
            name = line[:-9]
            # We make a new heading with the current name under the data that
            # matches the current annotation
            current_annotation = dict()
            # We add keep the dictonary accessible under 'current_annotation',
            # so we can keep adding new data to it without having to keep track
            # of where it belongs
            annotations[ann][name] = current_annotation
            continue

        linedata = parse_line(line)

        # If we got no data, we reached the end of the annotated section
        if not linedata:
            current_annotation = None
            continue

        # Lets get the name of the data
        name = linedata.pop('name')
        # If we are in an annotated section, we add it to the current annotation
        if current_annotation is not None:
            current_annotation[name] = linedata
        # Otherwise, we add the newfound annotation to the dictionary in case
        # we find a corresponding section later on.
        # The data from the current line we add directly to the output data
        else:
            annotation = linedata.pop('annotation')
            annotations[annotation] = linedata
            data[name] = linedata

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

    # Parsing the values
    values = values.strip().split()
    # Are there any values
    if not values:
        return data

    # If there is a single value
    if len(values) == 1:
        value = values.pop()
        # Is the value a percentage
        if value.endswith('%'):
            data['percentage'] = float(value[:-1])
        # Otherwise, it is an integer
        else:
            data['count'] = int(value)
    elif len(values) == 2:
        # If there are multiple values, they are in the format
        # count (percentage%)
        count = values[0]
        percentage = values[1]

        # The percentage should be in the format: (12.34% or 100% or 0%)
        assert re.fullmatch(r'[(]\d+\.?\d*%[)]', percentage)

        # Add the percentage and the count to the data
        data['percentage'] = float(percentage[1:-2])
        data['count'] = int(count)

    return data
