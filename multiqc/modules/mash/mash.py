#!/usr/bin/env python

""" MultiQC submodule to parse output from Samtools stats """

import logging
from multiqc.modules.base_module import BaseMultiqcModule
from collections import OrderedDict
from multiqc import config
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """ mash dist module """ # mash screen would be similar

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Mash', anchor='mash',
        href="https://mash.readthedocs.io/en/latest/index.html",
        info="Fast genome and metagenome distance estimation using MinHash.")

        # Find and load any Mash reports
        self.mashorganisms_data = dict()
        self.mashorganisms_keys = list()
        for myfile in self.find_log_files('mash'):
            self.mashorganisms_data.update({myfile['s_name'] : self.parse_mash(myfile)})

        if len(self.mashorganisms_data) == 0:
            raise UserWarning
        # Filter to strip out ignored sample names
        self.mashorganisms_data = self.ignore_samples(self.mashorganisms_data)
        log.info("Found {} logs".format(len(self.mashorganisms_data)))
        self.write_data_file(self.mashorganisms_data, 'multiqc_mash')

        # Plot the data
        self.add_section( plot = self.mash_plot() )

    def parse_mash(self, myfile):
        parsed_data = dict()
        for line in myfile['f'].splitlines():
            split_line = line.split("\t")
            if split_line[3] == '0' :
                # split_line[0] contains the full sequence name in a horrendous format, now to get Genus and species. It's not perfect, but it works 99.99% of the time
                split_split_line = split_line[0].split("-")[7]
                if split_split_line.startswith('_'):
                    split_split_line = split_split_line[1:]
                if '.' in split_split_line:
                    split_split_line = split_split_line.split('.')[0]
                final_split = " ".join(split_split_line.split("_")[:2])
                # get the keys
                if final_split not in self.mashorganisms_keys:
                    self.mashorganisms_keys.append(final_split)
                # count each occurent of that organism for each mash result
                if final_split not in parsed_data:
                    parsed_data[final_split] = 1
                else:
                    parsed_data[final_split] = parsed_data[final_split] + 1
        return(parsed_data)

    def mash_plot(self):
        config = {
            'id' : "Mash_Distances",
            'title': "Mash: Distances",
            'xlab': "Organism",
            'ylab': "Sample"
        }
        return bargraph.plot(self.mashorganisms_data, self.mashorganisms_keys, config)
