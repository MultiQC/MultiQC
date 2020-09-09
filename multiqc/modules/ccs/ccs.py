#!/usr/bin/env python

""" MultiQC module to parse output from CCS """

import logging
import re

from collections import OrderedDict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialse the parent object
        super(MultiqcModule, self).__init__(name='ccs', anchor='ccs',
                href='https://github.com/PacificBiosciences/ccs',
                info='Generate highly accurate single-molecule consensus '
                     'reads')

        # To store the mod data
        self.mod_data = dict()
        plotdata = dict()
        self.parse_log_files()
        self.write_data_files()
        self.add_general_stats()
        self.add_sections()

    def parse_log_files(self):
        for f in self.find_log_files('ccs', filehandles=True):
            data = parse_PacBio_log(f['f'])
            filename = f['s_name']
            self.mod_data[filename] = data
            self.add_data_source(f)

        # If we found no data
        if not self.mod_data:
            raise UserWarning

    def write_data_files(self):
        for filename in self.mod_data:
            self.write_data_file(self.mod_data[filename], filename)

    def add_general_stats(self):
        """ Add results overview to the general stats table """
        general_keys = ['ZMWs input', 'ZMWs generating CCS', 'ZMWs filtered']

        # Gather the data
        general_stats = dict()
        for filename in self.mod_data:
            data = self.mod_data[filename]
            file_data = OrderedDict()
            for key in general_keys:
                file_data[key] = data[key]['count']
            general_stats[filename] = file_data

        # Determine the formatting
        headers = OrderedDict()
        for i,field in enumerate(general_keys):
            headers[field] = {'placement': i}
        self.general_stats_addcols(general_stats, {key: {} for key in general_keys})

    def add_sections(self):
        plot_data = dict()
        # First we gather all the data we need
        for filename in self.mod_data:
            filter_reasons = self.filter_and_pass(self.mod_data[filename])
            plot_data[filename] = filter_reasons


        # Gather all the filter reasons we have found in the parsed reports
        reasons = set()
        for filename in plot_data:
            for reason in plot_data[filename]:
                reasons.add(reason)

        # Put them in a sorted list
        pass_ = 'ZMWs generating CCS'
        reasons = sorted(list(reasons))
        # Put ZMWS generating CCS at the front
        reasons.remove(pass_)
        reasons.insert(0,pass_)

        # Now we add the formatting we want
        # We want the passed ZMWs to be green
        green = '#5cb85c'
        formatting = OrderedDict()
        for reason in reasons:
            d = dict()
            d['name'] = reason
            # We want the PASS count to be green
            if reason == pass_:
                d['color'] = green
            formatting[reason] = d

        self.add_section (
                name='ZMWs filtered and passed',
                anchor='ccs-filter',
                description=(
                    'The number of reads that failed or passed all '
                    'filters'
                ),
                helptext=(
                    'The number of reads that passed all filters is shown '
                    'as **ZMWs generating CCS**. All other categories that '
                    'are shown in the graph represent the number of reads '
                    'that were dropped for the specified reason.'
                ),
                plot=bargraph.plot(plot_data, formatting)
        )

    def filter_and_pass(self, data):
        reasons = dict()

        # Add why ZMWs were filtered
        filter_data = data['ZMWs filtered']['Exclusive ZMW counts']
        for reason in filter_data:
            reasons[reason] = filter_data[reason]['count']
        # Add the number of ZMWs that passed
        reasons['ZMWs generating CCS'] = data['ZMWs generating CCS']['count']

        return reasons


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
        # If so, we will only add an empty dictionary with the correct name
        # These field are of the format "something something for (A):"
        section_header_pattern = ' for [(][A-Z][)][:]$'
        if re.search(section_header_pattern, line):
            linedata = parse_line(line)
            ann = linedata['annotation']
            # Cut off the ' for (B):' part
            name = line[:-9]
            # We make a new heading with the current name under the data that
            # matches the current annotation
            current_annotation = dict()
            # We keep the dictonary accessible under 'current_annotation',
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
        # If we are in an annotated section, we add the data to the current
        # annotation
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
    # Are there are no values we are done
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
