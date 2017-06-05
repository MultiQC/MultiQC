#!/usr/bin/env python

""" MultiQC module to parse results from jellyfish  """

from __future__ import print_function

from collections import OrderedDict
import logging
from multiqc import config
from multiqc.plots import linegraph, bargraph
from multiqc.modules.base_module import BaseMultiqcModule



# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='jellyfish', anchor='jellyfish',
        href="http://www.cbcb.umd.edu/software/jellyfish/",
        info="JELLYFISH is a tool for fast, memory-efficient counting of k-mers in DNA.")

        self.jellyfish_data_all = dict()
        for f in self.find_log_files(config.sp['jellyfish'], filehandles=True):
            self.parse_jellyfish_data(f)
        
        #clean histogram to prevent long flat tail
        self.clean_jellyfish_data()
        
        if len(self.jellyfish_data) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning
            
        log.info("Found {} reports".format(len(self.jellyfish_data)))
        
        self.intro += """Jellyfish k-mer plots for library complexity. 
        A possible way to assess the complexity of a
        library even in absence of a reference sequence is to
        look at the kmer profile of the reads. The idea is to
        count all the kmers (i.e., sequence of length k) that occur
        in the reads. In this way it is possible to know how many
        kmers occur 1,2,..., N times and represent this as a
        plot. This plot tell us for each x, how many k-mers
        (y-axis) are present in the dataset in exactly x-copies.
        In an ideal world (no errors in sequencing, no bias, no
        repeated regions) this plot should be as close as
        possible to a gaussian distribution. In reality we will
        always see a peak for x=1 (i.e., the errors) and another
        peak close to the expected coverage. If the genome is
        highly heterozygous a second peak at half of the coverage
        can be expected."""
        
        
        self.sections = list()
        self.sections.append({
            'name': 'Jellyfish plot for k-mers with coverage between 0 and {}'.format(self.jellyfish_max_x),
            'anchor': 'jellyfish_kmer_plot',
            'content': self.frequencies_plot()
        })

    def clean_jellyfish_data(self):
        """ Avoid to plot long flat tails by loosing 0,5% of the tail """
        max_x = 0
        #find where the max is
        for s_name, data in self.jellyfish_data_all.items():
            max_key = max(data, key=data.get)
            if max_key < 100:
                max_key = 200 # the maximum is below 100, we display anyway up to 200
            else:
                max_key = 2*max_key #in this case the area plotted is a function of the maximum x
            max_x = max(max_x, max_key)
        self.jellyfish_max_x = max_x
        #create the container I am going to plot
        self.jellyfish_data = dict()
        for s_name, d in self.jellyfish_data_all.items():
            self.jellyfish_data[s_name] = dict([(i,d[i]) for i in range(max_x+1)])


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
        if len(histogram) > 0:
            if f['s_name'] in self.jellyfish_data_all:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            self.add_data_source(f)
            self.jellyfish_data_all[f['s_name']] = histogram



    def frequencies_plot(self, xmin=0):
        """ Generate the qualities plot """
        pconfig = {
            'id': 'Jellyfish_kmer_plot',
            'title': 'Jellyfish: K-mer plot',
            'ylab': 'Counts',
            'xlab': 'k-mer frequency',
            'xDecimals': False,
            'xmin': xmin
        }
        return linegraph.plot(self.jellyfish_data, pconfig)



