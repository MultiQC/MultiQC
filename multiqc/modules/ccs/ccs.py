#!/usr/bin/env python

""" MultiQC module to parse output from CCS """

import json
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
        super(MultiqcModule, self).__init__(name='CCS', anchor='ccs',
                href='https://github.com/PacificBiosciences/ccs',
                info=(
                    'is used to generate highly accurate single-molecule '
                    'consensus reads'
                )
        )

        # To store the mod data
        self.mod_data = dict()
        self.parse_v4_log_files()
        self.parse_v5_log_files()
        self.write_data_files()
        self.add_sections()

    def parse_v4_log_files(self):
        for f in self.find_log_files('ccs/v4', filehandles=True):
            data = parse_PacBio_log(f['f'])
            v5_data = self.convert_to_v5(data)
            filename = f['s_name']
            self.mod_data[filename] = v5_data
            self.add_data_source(f)

    def parse_v5_log_files(self):
        for f in self.find_log_files('ccs/v5', filehandles=True):
            v5_data = json.load(f['f'])
            filename = f['s_name']
            self.mod_data[filename] = v5_data
            self.add_data_source(f)

        # If we found no data
        if not self.mod_data:
            raise UserWarning

    def write_data_files(self):
        self.write_data_file(self.mod_data, 'multiqc_ccs_report')

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
        pass_ = 'ZMWs pass filters'
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

        # Plot configuration
        config = {
                'id': 'ccs-filter-graph',
                'title': 'CCS: ZMW results',
                'ylab': 'Number of ZMWs',
                'xlab': 'CCS report file'
        }
        self.add_section (
                name='ZMWs filtered and passed',
                anchor='ccs-filter',
                description=(
                    'The number of ZMWs that failed or passed all '
                    'filters'
                ),
                helptext=(
                    'The number of ZMWs that passed all filters is shown '
                    'as **ZMWs generating CCS**. All other categories that '
                    'are shown in the graph represent the number of ZMWs '
                    'that were dropped for the specified reason.'
                ),
                plot=bargraph.plot(plot_data, formatting, config)
        )

    def filter_and_pass(self, data):
        """ Gather the reasons why ZMWs were failed or passed """
        reasons = dict()

        # We only have to use the attributes
        attributes = data['attributes']

        # Add filtere reasons (id starts with filtered_) and total passed
        for entry in attributes:
            if entry['id'].startswith('filtered') or entry['id'] == 'zmw_passed_yield':
                reasons[entry['name']] = entry['value']

        return reasons

    def convert_to_v5(self, data):
        """ Convert the v4 format to the new CCS v5 json format """
        # Initialise the v5 format dictionary
        v5 = dict()
        v5['id'] = 'ccs_processing'
        attributes = list()
        v5['attributes'] = attributes

        # Update names for top level entries, they have been changed in v5
        data['ZMWs pass filters'] = data['ZMWs generating CCS']
        data['ZMWs fail filters'] = data['ZMWs filtered']

        # Add the top level values to the attributes
        for name in ('ZMWs input', 'ZMWs pass filters', 'ZMWs fail filters'):
            count = data[name]['count']
            attributes.append(self.make_v5_entry(name, count))

        # Add the reasons for filtering to the attributes
        for reason in data['ZMWs filtered']['Exclusive ZMW counts']:
            count = data['ZMWs filtered']['Exclusive ZMW counts'][reason]['count']
            attributes.append(self.make_v5_entry(reason, count))

        return v5

    def make_v5_entry(self, name, count):
        """ Make a v5 output entry based on name and count """
        # Dictionary to map the report names to the v5 id annotation
        name_to_id = {
                'ZMWs input': 'zmw_input',
                'ZMWs pass filters': 'zmw_passed_yield',
                'ZMWs fail filters': 'zmw_filtered_yield',
                'Below SNR threshold': 'filtered_poor_snr',
                'Median length filter': 'filtered_median_length_filter',
                'Lacking full passes': 'filtered_too_few_passes',
                'Heteroduplex insertions': 'filtered_heteroduplexes',
                'Coverage drops': 'filtered_coverage_drops',
                'Insufficient draft cov': 'filtered_insufficient_draft_cov',
                'Draft too different': 'filtered_draft_too_different',
                'Draft generation error': 'filtered_draft_failure',
                'Draft above --max-length': 'filtered_draft_too_long',
                'Draft below --min-length': 'filtered_draft_too_short',
                'Reads failed polishing': 'filtered_read_failed_polish',
                'Empty coverage windows': 'filtered_empty_window_during_polishing',
                'CCS did not converge': 'filtered_non_convergent',
                'CCS below minimum RQ': 'filtered_below_rq',
                'Unknown error': 'filtered_unknown_error'
        }

        return {
                'name': name,
                'id': name_to_id[name],
                'value': count
                }


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

        # The percentage should be in the format:
        # (12.34%) or (100%) or (0%) or (-nan%)
        # So we can remove the bracets first
        percentage = percentage[1:-1]
        # Then we make sure it is one of the supported patterns
        assert re.fullmatch(r'(\d+\.\d+%|\d+%|-nan%)', percentage)

        # If the percentage is this weird nan, set it to 0
        if percentage == '-nan%':
            data['percentage'] = 0.0
        else:  # Otherwise, we cut of the % sign and convert to float
            data['percentage'] = float(percentage[:-1])
        # Add the count as an integer
        data['count'] = int(count)

    return data
