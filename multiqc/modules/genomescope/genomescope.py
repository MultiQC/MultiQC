#!/usr/bin/env python

""" MultiQC module to parse output from genomescope """

from __future__ import print_function
from collections import OrderedDict
import logging
import re
from multiqc.plots import table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """Genomescope module"""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="GENOMESCOPE",
            anchor="genomescope",
            href="https://github.com/schatzlab/genomescope",
            info="Fast genome analysis from unassembled short reads.",
            doi="https://doi.org/10.1093/bioinformatics/btx153",
        )


        # Find and load any MERQURY reports
        self.summary_data = dict()
        for f in self.find_log_files("genomescope/summary"):
            self.parse_summary_log(f)

        log.info("Found {} summary reports".format(len(self.summary_data)))

        # Write parsed report data to a file
        self.write_data_file(self.summary_data, "multiqc_genomescope_summary")
        
        
        #self.merqury_general_stats_table()
        self.add_section(name="Summary", anchor="genomescope-summary", plot=self.summary_table())
    
    def parse_summary_log(self, f):
        self.add_data_source(f,f["s_name"])
        s_name=f["s_name"].replace("_Summary.txt","")
        block=-1
        i=0
        content=[]
        for l in f["f"].splitlines():
            if l.startswith("property"): block=i
            if block!=-1: content.append(l)
            if l.startswith("GenomeScope version "): version=l.replace("GenomeScope version ","")
            if l.startswith("k = "): k=l.split(" = ")[1]
            if l.startswith("p = "): p=l.split(" = ")[1]
            i=i+1
        legend=content[0].split()
        self.summary_data[s_name]= dict()
        self.legends=[]
        for l in content[1:]:
            items= re.split(r'\s{2,}', l.strip())
            this_name=items[0]
            this_val=items[2]
            self.summary_data[s_name][this_name]= this_val
            self.legends.append(this_name)
        
    def summary_table(self):
        """Take the parsed stats from the QUAST report and add some to the
        General Statistics table at the top of the report"""
        headers = OrderedDict()
        for legend in self.legends:
            headers[legend] = {
                "title": legend,
                "description": legend,
                "min": 0,
                #"suffix": self.contig_length_suffix,
                "scale": "RdYlGn",
                "format": "{:,.2f}",
            }

        config = {
            "id": "summary_table",
            "namespace": "GENOMESCOPE",
            "min": 0,
        }
        return table.plot(self.summary_data,headers,config)
                

