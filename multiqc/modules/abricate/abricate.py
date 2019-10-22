#!/usr/bin/env python

""" MultiQC module to parse output from ABRicate """

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import heatmap
import logging
import re

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='ABRicate', anchor='abricate',
        href="https://github.com/tseemann/abricate",
        info="Mass screening of contigs for antimicrobial resistance or virulence genes.")

        self.abricate_data  = dict()
        self.abricate_xcats = dict()
        self.abricate_ycats = dict()
        # Find all files for mymod
        for myfile in self.find_log_files('abricate'):
            db = myfile['s_name']
            db = db.replace('.', '_') # just in case the user puts a period in the name, HTML hates that
            self.getdata(myfile, db)
            self.add_section(
                name = db,
                anchor = 'abricate-' + db,
                description = 'This heatmap shows the score for each sample with the genes supplied from ' + db,
                plot = self.abricate_heatmap_plot(db) )

        log.info("Found {} logs".format(len(self.abricate_data)))

    def getdata(self, myfile, db):
        self.abricate_ycats[db] = []
        self.abricate_data[db] = []
        for line in myfile['f'].splitlines():
            line = line.replace('\t.', '\t0.00')
            if not line.split("\t")[0] == "#FILE":
                # gets the sample name
                self.abricate_ycats[db].append(line.split("\t")[0])
                # gets the numbers
                # sometimes organisms are found to carry multiple copies, so those elements need to be converted to a single number
                if ';' in line:
                    double_carrier_line = line.split("\t")[2:]
                    for status in double_carrier_line:
                        if ';' in status:
                            # get the max number of multiple hits
                            double_carrier_line[double_carrier_line.index(status)]=max(status.split(';'), key=float)
                    self.abricate_data[db].append(double_carrier_line)
                else:
                    self.abricate_data[db].append(line.split("\t")[2:])
            else:
                # gets the gene names
                self.abricate_xcats[db] = line.split("\t")[2:]

    def abricate_heatmap_plot(self, db):
        config = {
            'id' : "abricate_" + db,
            'title': "ABRicate: " + db,
            'square': False,
            'colstops': [ [0, '#FFFFFF'], [0.6, '#ffffe5'], [0.7, '#d9f0a3'], [0.95, '#004529'], [1, '#000000'], ]
        }
        return heatmap.plot(self.abricate_data[db], self.abricate_xcats[db], self.abricate_ycats[db], config)
