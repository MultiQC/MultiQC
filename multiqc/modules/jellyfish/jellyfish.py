#!/usr/bin/env python

""" MultiQC module to parse results from jellyfish  """

from __future__ import print_function

from collections import OrderedDict
import logging
from multiqc import config, BaseMultiqcModule, plots

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
        plot. This plot tell us for each x, how many 25-mers
        (y-axis) are present in the dataset in exactly x-copies.
        In an ideal world (no errors in sequencing, no bias, no
        repeating regions) this plot should be as close as
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
        if self.detect_high_unique_kmers():
            self.sections.append({
                'name': 'Jellyfish plot for k-mers with coverage between 3 and {}'.format(self.jellyfish_max_x),
                'anchor': 'jellyfish_kmer_plot_no_low_freq',
                'content': 'It has been detected that the number of unique k- mers is unexpectedly high. This is not a good sign (it means a lot of unique and most likely erroneus k-mers). This plot removes the k-mers occuring in single and double copy to allow a better visualisation of the plot.' + self.frequencies_plot(10)
            })

    def clean_jellyfish_data(self):
        """ Avoid to plot long flat tails by loosing 0,5% of the tail """
        max_x = 0
        total_bases_by_sample = dict()
        for s_name, d in self.jellyfish_data_all.items():
            total_bases_by_sample[s_name] = sum(d.values())
            cumulative = 0
            for count in sorted(d.keys(), reverse=True):
                cumulative += d[count]
                if cumulative / float(total_bases_by_sample[s_name]) > 0.05:
                    max_x = max(max_x, count)
                    break
        self.jellyfish_max_x = max_x
        #create the container I am going to plot
        self.jellyfish_data = dict()
        for s_name, d in self.jellyfish_data_all.items():
            self.jellyfish_data[s_name] = dict([(i,d[i]) for i in range(max_x+1)])

    def detect_high_unique_kmers(self):
        """ sometime there are a lot of unique kmers. This indicates a bad library but some time can be a property of the data it self  """
        return True

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
        

        if len(histogram) > 0:
            if f['s_name'] in self.jellyfish_data_all:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            self.add_data_source(f)
            self.jellyfish_data_all[f['s_name']] = histogram



    def frequencies_plot(self, xmin=0):
        """ Generate the qualities plot """
        
        pconfig = {
            'smooth_points': 200,
            'id': 'Jellyfish_kmer_plot',
            'title': 'Jellyfish: K-mer plot',
            'ylab': 'Count',
            'xlab': 'Frequency',
            'xDecimals': False,
            'ymin': 0,
            'xmin': xmin
        }

        return plots.linegraph.plot(self.jellyfish_data, pconfig)



