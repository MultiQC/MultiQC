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

        self.jellyfish_data = dict()
        for f in self.find_log_files(config.sp['jellyfish'], filehandles=True):
            self.parse_jellyfish_data(f)
        
        if len(self.jellyfish_data) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning
            
        log.info("Found {} reports".format(len(self.jellyfish_data)))
         
        self.intro += 'My amazing module output'
        self.intro += self.frequencies_plot()
        

    def parse_jellyfish_data(self, f):
        """ Go through the hist file and memorise it """
        histogram = {}
        for line in f['f']:
            line = line.rstrip('\n')
            occurence = int(line.split(" ")[0])
            count = int(line.split(" ")[1])
            histogram[occurence] = occurence*count

        #### this nees to go in init or in another funciotn as phil did
        #### remove the tail by using this code (modify it)
        #### if position 1 is 5 times higher than max of pos[2...max] then plot two plots and a warning 
        if len(histogram) > 0:
        # Chew back on histogram to prevent long flat tail
        # (find a sensible max x - lose 1% of longest tail)
            max_x = 0
            import pdb
            pdb.set_trace()
            sum_of_occurces = sum(historgra.values())
            cumulative = 0
            for count in sorted(d.keys(), reverse=True):
                cumulative += d[count]
                if cumulative / total_bases_by_sample[s_name] > 0.01:
                    max_x = max(max_x, count)
                    break


        if len(histogram) > 0:
            if f['s_name'] in self.jellyfish_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(f['s_name']))
            self.add_data_source(f)
            self.jellyfish_data[f['s_name']] = histogram



    def frequencies_plot(self):
        """ Generate the qualities plot """
        
        pconfig = {
            'smooth_points': 200,
            'id': 'Jellyfish_kmer_plot',
            'title': 'Jellyfish: K-mer plot',
            'ylab': 'Count',
            'xlab': 'Frequency',
            'xDecimals': False,
            'ymin': 0
        }

        return plots.linegraph.plot(self.jellyfish_data, pconfig)



