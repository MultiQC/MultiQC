#!/usr/bin/env python

""" MultiQC module to parse output from SortMeRNA """

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import re

from multiqc import config
from multiqc.plots import bargraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='SortMeRNA', anchor='sortmerna',
        href='http://bioinfo.lifl.fr/RNA/sortmerna/',
        info="is a program for filtering, mapping and OTU-picking NGS reads in metatranscriptomic and metagenomic data.")

        # Parse logs
        self.sortmerna = dict()
        for f in self.find_log_files('sortmerna', filehandles=True):
            self.parse_sortmerna(f)

        # Filter to strip out ignored sample names
        self.sortmerna = self.ignore_samples(self.sortmerna)

        if len(self.sortmerna) == 0:
            log.debug("Could not find any SortMeRNA data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} logs".format(len(self.sortmerna)))
        self.write_data_file(self.sortmerna, 'multiqc_sortmerna')

        # Get custom table header, default to 'rRNA'
        tab_header = getattr(config, 'sortmerna', {}).get('seqname', '% rRNA')

        # Add rRNA rate to the general stats table
        headers = OrderedDict()
        headers['rRNA_pct'] = {
            'title': tab_header,
            'description': 'Percentage of reads matched to a SortMeRNA database',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'OrRd'
        }
        self.general_stats_addcols(self.sortmerna, headers)

        # Make barplot
        self.sortmerna_detailed_barplot()

    def parse_sortmerna(self, f):
        s_name = None
        post_results_start = False
        post_database_start = False
        db_number = 0
        err = False

        for l in f["f"]:
            if "Reads file" in l:
                s_name = self.clean_s_name( l.split('=')[-1], f['root'] )
                self.sortmerna[s_name] = dict()
            if "Results:" in l and not post_results_start:
                post_results_start = True
            if not post_results_start:
                continue
            if post_results_start and not post_database_start:
                if "Total reads =" in l:
                    m = re.search("\d+",l)
                    if m:
                        self.sortmerna[s_name]["total"] = int(m.group())
                    else:
                        err = True
                elif "Total reads passing" in l:
                    m = re.search("\d+",l)
                    if m:
                        self.sortmerna[s_name]["rRNA"] = int(m.group())
                        self.sortmerna[s_name]["rRNA_pct"] = float(self.sortmerna[s_name]["rRNA"]) / float(self.sortmerna[s_name]["total"]) * 100
                    else:
                        err = True
                elif "Total reads failing" in l:
                    m = re.search("\d+",l)
                    if m:
                        self.sortmerna[s_name]["non_rRNA"] = int(m.group())
                        self.sortmerna[s_name]["non_rRNA_pct"] = float(self.sortmerna[s_name]["non_rRNA"]) / float(self.sortmerna[s_name]["total"]) * 100
                    else:
                        err = True
            if post_database_start:
                if not l.strip():
                    break
                db_number = db_number + 1
                m = re.search("    .*\t", l)
                if m:
                    db = m.group().strip()
                    db = os.path.splitext(os.path.basename(db))[0]
                    pct = float(re.search("\d+\.\d+%", l).group().replace("%",""))
                    count = int(self.sortmerna[s_name]["total"]) * (pct / 100)
                    self.sortmerna[s_name][db + "_pct"] = pct
                    self.sortmerna[s_name][db + "_count"] = count
                else:
                    err = True
            if "By database:" in l and not post_database_start:
                post_database_start = True
        if err:
            log.warning("Error parsing data in: " + s_name)
            self.sortmerna.pop(s_name, 'None')
        s_name = None

    def sortmerna_detailed_barplot (self):
        """ Make the HighCharts HTML to plot the sortmerna rates """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        metrics = set()
        for sample in self.sortmerna:
            for key in self.sortmerna[sample]:
                if not key in ["total", "rRNA", "non_rRNA"] and not "_pct" in key:
                    metrics.add(key)

        for key in metrics:
            keys[key] = { 'name': key.replace("_count","") }

        # Config for the plot
        pconfig = {
            'id': 'SortMeRNA detailed',
            'title': 'SortMeRNA hits',
            'ylab': 'Reads'
        }

        self.add_section( plot = bargraph.plot(self.sortmerna, keys, pconfig) )
