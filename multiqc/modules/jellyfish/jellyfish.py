#!/usr/bin/env python

""" MultiQC module to parse results from jellyfish  """

from __future__ import print_function

import logging
from multiqc import config
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Jellyfish', anchor='jellyfish',
        href="http://www.cbcb.umd.edu/software/jellyfish/",
        info="is a tool for fast, memory-efficient counting of k-mers in DNA.")

        self.jellyfish_data  = dict()
        self.jellyfish_max_x = 0
        for f in self.find_log_files('jellyfish', filehandles=True):
            self.parse_jellyfish_data(f)

        if self.jellyfish_max_x < 100:
            self.jellyfish_max_x = 200 # the maximum is below 100, we display anyway up to 200
        else:
            self.jellyfish_max_x = 2*self.max_key #in this case the area plotted is a function of the maximum x

        if len(self.jellyfish_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.jellyfish_data)))

        self.frequencies_plot(xmax=self.jellyfish_max_x)


    def parse_jellyfish_data(self, f):
        """ Go through the hist file and memorise it """
        histogram = {}
        occurence = 0
        for line in f['f']:
            line = line.rstrip('\n')
            occurence = int(line.split(" ")[0])
            count = int(line.split(" ")[1])
            histogram[occurence] = occurence*count
        #delete last occurnece as it is the sum of all kmer occuring more often than it.
        del histogram[occurence]
        #sanity check
        self.max_key  = max(histogram, key=histogram.get)
        self.jellyfish_max_x = max(self.jellyfish_max_x, self.max_key)
        if len(histogram) > 0:
            if f['s_name'] in self.jellyfish_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            self.add_data_source(f)
            self.jellyfish_data[f['s_name']] = histogram


    def frequencies_plot(self, xmin=0, xmax=200):
        """ Generate the qualities plot """

        helptext = '''
            A possible way to assess the complexity of a library even in
            absence of a reference sequence is to look at the kmer profile of the reads.
            The idea is to count all the kmers (_i.e._, sequence of length `k`) that occur
            in the reads. In this way it is possible to know how many kmers occur
            `1,2,.., N` times and represent this as a plot.
            This plot tell us for each x, how many k-mers (y-axis) are present in the
            dataset in exactly x-copies.

            In an ideal world (no errors in sequencing, no bias, no  repeated regions)
            this plot should be as close as  possible to a gaussian distribution.
            In reality we will always see a peak for `x=1` (_i.e._, the errors)
            and another peak close to the expected coverage. If the genome is highly
            heterozygous a second peak at half of the coverage can be expected.'''

        pconfig = {
            'id': 'Jellyfish_kmer_plot',
            'title': 'Jellyfish: K-mer plot',
            'ylab': 'Counts',
            'xlab': 'k-mer frequency',
            'xDecimals': False,
            'xmin': xmin,
            'xmax': xmax
        }

        self.add_section(
            anchor = 'jellyfish_kmer_plot',
            description = 'The K-mer plot lets you estimate library complexity and coverage from k-mer content.',
            helptext = helptext,
            plot = linegraph.plot(self.jellyfish_data, pconfig)
        )
