#!/usr/bin/env python

""" MultiQC module to parse output from WhatsHap """

from collections import defaultdict
import logging

from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    WhatsHap module class, parses WhatsHap output.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
                name='WhatsHap', anchor='whatshap',
                href='https://whatshap.readthedocs.io/',
                info=(
                    'WhatsHap is a software for phasing genomic variants '
                    'using sequencing reads.'
                )
        )

        self.whatshap_stats = dict()

        for f in self.find_log_files('whatshap/stats', filehandles=True):
            sample, data = self.parse_whatshap_stats(f['f'])
            self.whatshap_stats[sample] = data

        # Filter to strip out ignored sample names
        self.whatshap_stats = self.ignore_samples(self.whatshap_stats)

        if len(self.whatshap_stats) == 0:
            raise UserWarning

        # Write parsed report data to a file
        self.write_data_file(self.whatshap_stats, 'multiqc_whatshap_stats')


    def parse_whatshap_stats(self, file_content):
        """ Parse WhatsHap stats file """

        def process_data(data):
            """
            Process the data to convert it to the proper format

            In practice, this only involves the following changes:
                - Convert all numeric values to integers
                - Replace 'nan' values to python None
            """
            for key, value in data.items():
                # Try to make the value an int
                try:
                    data[key] = int(data[key])
                except ValueError:
                    pass

                # Replace 'nan' with None
                if value == 'nan':
                    data[key] = None
                
        # Get the header and remove the "#" character from the first field
        header = next(file_content).strip().split()
        header[0] = header[0][1:]

        # The results from this log file, a dictionary for every sample
        results = defaultdict(dict)
        # Parse the lines that have the data
        for line in file_content:
            spline = line.strip().split()
            data = {field:value for field, value in zip(header, spline)}
            process_data(data)
            # Remove the sample and chromsome from the data
            sample = data.pop('sample')
            chromosome = data.pop('chromosome')
            # Insert the current line under chromosome
            results[chromosome] = data
        return sample, results
