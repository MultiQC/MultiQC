#!/usr/bin/env python

""" MultiQC module to parse output from Lima """

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
        super(MultiqcModule, self).__init__(name='Lima', anchor='lima',
              href='https://github.com/PacificBiosciences/barcoding',
              info= ' is used to demultiplex single-molecule sequencing reads.'
        )

        # To store the summary data
        self.lima_summary = dict()
        self.lima_counts = dict()
        self.parse_summary_files()
        self.parse_counts_files()
        self.write_data_files()
        # Add the counts to the general statistics
        self.general_stats_addcols(self.lima_counts)

    def parse_summary_files(self):
        for f in self.find_log_files('lima/summary', filehandles=True):
            data = parse_PacBio_log(f['f'])
            filename = f['s_name']
            self.lima_summary[filename] = data
            self.add_data_source(f)

    def parse_counts_files(self):
        for f in self.find_log_files('lima/counts', filehandles=True):
            data = self.parse_lima_counts(f['f'], f['root'])
            # Check for duplicate samples
            for sample in data:
                if sample in self.lima_counts:
                    msg = f'Duplicate sample name found! Overwriting: {sample}'
                    log.debug(msg)
            # After warning the user, overwrite the samples
            self.lima_counts.update(data)
            self.add_data_source(f)
        # Remove samples that were specified using --ignore-samples
        self.lima_counts = self.ignore_samples(self.lima_counts)

    def write_data_files(self):
        if not self.lima_summary and not self.lima_counts:
            raise UserWarning
        if self.lima_summary:
            self.write_data_file(self.lima_summary, 'multiqc_lima_summary')
        if self.lima_counts:
            self.write_data_file(self.lima_counts, 'multiqc_lima_counts')

    def parse_lima_counts(self, file_content, root):
        """ Parse lima counts file """

        # The first line is the header
        header = next(file_content).strip().split()

        # A dictionary to store the results
        lima_counts = dict()
        for line in file_content:
            spline = line.strip().split()
            data = {field: value for field, value in zip(header, spline)}

            sample = self.clean_s_name(data['IdxFirstNamed'], root)
            counts = data['Counts']
            mean_score = data['MeanScore']
            lima_counts[sample] = {
                    'Counts': counts,
                    'MeanScore': mean_score
            }

        return lima_counts


def parse_PacBio_log(file_content):
    """ Parse PacBio log file """
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
    """ Parse a line from the Lima log file """
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
