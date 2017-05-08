#!/usr/bin/env python

""" MultiQC module to parse output from Samblaster """

from __future__ import print_function
import os
from collections import OrderedDict
import logging
import re
from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """ Samblaster """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Samblaster', anchor='samblaster',
                                            href="https://github.com/GregoryFaust/samblaster",
                                            info="is a tool to mark duplicates and extract discordant and split reads from sam files.")

        self.samblaster_data = dict()
        for f in self.find_log_files('samblaster', filehandles=True):
            self.parse_samblaster(f)

        # Filter to strip out ignored sample names
        self.samblaster_data = self.ignore_samples(self.samblaster_data)

        if len(self.samblaster_data) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning

        headers = OrderedDict()
        headers['pct_dups'] = {
            'title': '% Dups',
            'description': 'Percent Duplication',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'OrRd'
        }

        self.general_stats_addcols(self.samblaster_data, headers)

        # Write parsed report data to a file
        self.write_data_file(self.samblaster_data, 'multiqc_samblaster')

        log.info("Found {} reports".format(len(self.samblaster_data)))

        self.add_barplot()

    def add_barplot(self):
        """ Generate the Samblaster bar plot. """
        cats = OrderedDict()
        cats['n_nondups'] = {'name': 'Non-duplicates'}
        cats['n_dups'] = {'name': 'Duplicates'}

        pconfig = {
            'id': 'samblaster_duplicates',
            'title': 'Number of duplicate reads',
        }
        self.add_section( plot = bargraph.plot(self.samblaster_data, cats, pconfig) )

    def parse_samblaster(self, f):
        """ Go through log file looking for samblaster output.
        If the
        Grab the name from the RG tag of the preceding bwa command """
        dups_regex = "samblaster: (Removed|Marked) (\d+) of (\d+) \((\d+.\d+)%\) read ids as duplicates"
        input_file_regex = "samblaster: Opening (\S+) for read."
        rgtag_name_regex = "\\\\tID:(\S*?)\\\\t"
        data = {}
        s_name = None
        fh = f['f']
        for l in fh:
            # try to find name from RG-tag. If bwa mem is used upstream samblaster with pipes, then the bwa mem command
            # including the read group will be written in the log
            match = re.search(rgtag_name_regex, l)
            if match:
                s_name = match.group(1)

            # try to find name from the input file name, if used
            match = re.search(input_file_regex, l)
            if match:
                basefn = os.path.basename(match.group(1))
                fname, ext = os.path.splitext(basefn)
                # if it's stdin, then try bwa RG-tag instead
                if fname != 'stdin':
                    s_name = fname

            match = re.search(dups_regex, l)
            if match:
                data['n_dups'] = int(match.group(2))
                data['n_tot'] = int(match.group(3))
                data['n_nondups'] = data['n_tot'] - data['n_dups']
                data['pct_dups'] = float(match.group(4))

        if s_name is not None:
            s_name = self.clean_s_name(s_name, f['root'])
            if s_name in self.samblaster_data:
                log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f['fn'], s_name))
            self.add_data_source(f, s_name)
            self.samblaster_data[s_name] = data
