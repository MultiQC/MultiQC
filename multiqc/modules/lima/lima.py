#!/usr/bin/env python

""" MultiQC module to parse output from Lima """

import logging
import re

from collections import OrderedDict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph
from multiqc.utils import config

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
        # Add a graph of all filtered ZMWs
        self.add_sections()
        # Add a plot for the data values in the counts file
        self.plot_counts_data()

    def parse_summary_files(self):
        for f in self.find_log_files('lima/summary', filehandles=True):
            data = parse_PacBio_log(f['f'])
            filename = f['s_name']
            self.lima_summary[filename] = data
            self.add_data_source(f)

    def parse_counts_files(self):
        for f in self.find_log_files('lima/counts', filehandles=True):
            data = self.parse_lima_counts(f['f'], f['root'])
            # Update sample names if --lima-barcodes was specified
            self.update_sample_names(data)
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

            first_barcode = data['IdxFirstNamed']
            second_barcode = data['IdxCombinedNamed']
            # The format barcode1--barcode2 is also used in
            # self.update_sample_names
            sample = self.clean_s_name(f'{first_barcode}--{second_barcode}', root)
            counts = data['Counts']
            mean_score = data['MeanScore']
            lima_counts[sample] = {
                    'Counts': int(counts),
                    'MeanScore': float(mean_score)
            }

        return lima_counts

    def update_sample_names(self, data):
        lima_barcodes = config.kwargs['lima_barcodes']
        # If --lima-barcodes wasn't specified
        if not lima_barcodes:
            return

        # Read the barcodes to sample mapping
        barcodes = self.read_lima_barcodes(lima_barcodes)

        # Copy the data over to the new dictionary by sample name
        for name in list(data.keys()):
            # If we have a mapping from the barcode to samplename
            if name in barcodes:
                # Get the new sample name
                sample_name = barcodes[name]
                # Update the sample name
                data[sample_name] = data.pop(name)

    def read_lima_barcodes(self, lima_barcodes):
        """
        Read the lima barcodes file, and return as a dictionary

        The keys will be of the format barcode1--barcode2, and the values will
        be the sample name
        """
        barcodes = dict()
        with open(lima_barcodes) as fin:
            for line in fin:
                sample, first_barcode, second_barcode = line.strip().split('\t')
                barcodes[f'{first_barcode}--{second_barcode}'] = sample
        return barcodes

    def add_sections(self):
        plot_data = dict()

        # First, we gather all filter results for each lima summary
        all_filters = set()
        for filename, data in self.lima_summary.items():
            for reason in self.filter_and_pass(data):
                all_filters.add(reason)

        plot_data = list()
        for filename, data in self.lima_summary.items():
            d = dict()
            for reason in all_filters:
                d[reason] = dict()
                # This is where the filter reasons are stored
                zmw_marginals = data['ZMWs below any threshold']['ZMW marginals']
                for filter_reason in zmw_marginals:
                    if filter_reason == reason:
                        d[reason][filename] = zmw_marginals[filter_reason]['count']
                        break
                # If we didn't find filter reason for this filename, set it to
                # zero
                else:
                    d[reason][filename] = 0
                # Unless 'reason' is actually the special case of the number of
                # unfiltered reads
                if reason == 'ZMWs above all thresholds':
                    d[reason][filename] = data['ZMWs above all thresholds']['count']
            plot_data.append(d)

        # Plot configuration
        config = {
                'id': 'lima-filter-graph',
                'title': 'Lima: ZMW results',
                'ylab': 'Number of ZMWs',
                'xlab': 'Lima summary file',
                'cpswitch': False,
                'logswitch': True,
                'stacking': None,
                'data_labels': [
                    {
                        'name': filename,
                        'cpswitch_counts_label': filename
                    } for filename in self.lima_summary
                ]
        }

        formatting = [filename for filename in self.lima_summary]

        self.add_section(
                name='ZMWs filtered by Lima',
                anchor='lima-filter',
                description=(
                    'The number of ZMWs that failed or passed all Lima '
                    'filters'
                ),
                helptext=(
                    'The number of ZMWs that passed all filters is shown as '
                    '**ZMWs above all thresholds**, all other categories that '
                    'are shown in the graph represent the number of ZMWs that '
                    'were dropped for the specified reason.'
                ),
                plot = bargraph.plot(plot_data, formatting, config)
        )

    def filter_and_pass(self, data):
        """ Get the reasons why each ZMW was filtered """
        reasons = dict()

        # Add why ZMWs were filtered
        filter_data = data['ZMWs below any threshold']['ZMW marginals']
        for reason in filter_data:
            reasons[reason] = filter_data[reason]['count']
        # Add the ZMWs that passed
        reasons['ZMWs above all thresholds'] = data['ZMWs above all thresholds']['count']

        return reasons

    def plot_counts_data(self):
        """ Plot the lima counts data """
        pdata = list()

        categories = ['Counts', 'MeanScore']
        description = ['Number of Reads', 'Mean Quality Score']

        configuration = {
            'id': 'multiqc_lima_counts',
            'title': 'Lima: Number of Reads',
            'anchor': 'multiqc_lima_counts',
            # Placeholder y-label, the actual values will be set below in
            # data_labels
            'ylab': '# Reads',
            'data_labels': [
                {
                    'name': category,
                    'cpswitch_counts_label': category,
                    'ylab': description
                } for category, description in zip(categories, description)
            ]
        }

        for category in categories:
            data = dict()
            for sample, entry in self.lima_counts.items():
                data[sample] = { category: entry[category] }
            pdata.append(data)

        self.add_section(
                name='Per sample count data',
                anchor='multiqc_lima_count',
                description=
                """
                    Per sample or per barcode statistics from Lima.

                    The **Counts** show the
                    number of reads for each sample or barcode pair.

                    The **MeanScore** shows the mean quality score of the reads for
                    each sample or barcode pair.
                """,
                helptext=
                """
                To display sample names instead of `barcode--barcode` pairs,
                please re-run MultiQC and specify which barcodes belong to
                which sample using the `--lima-barcodes` flag.
                """,

                plot=bargraph.plot(pdata, categories, configuration)),


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
