#!/usr/bin/env python

""" MultiQC module to parse output from SortMeRNA """

from __future__ import print_function
from collections import OrderedDict
import json
import logging
import os
import re

from multiqc.modules.base_module import config, BaseMultiqcModule
from multiqc import plots

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='SortMeRNA', anchor='sortmerna',
        href='http://bioinfo.lifl.fr/RNA/sortmerna/',
        info="is a program tool for filtering, mapping and OTU-picking NGS reads in metatranscriptomic and metagenomic data.")

        # Parse logs
        self.sortmerna = dict()
        for f in self.find_log_files(config.sp['sortmerna'], filehandles=True):
            self.parse_sortmerna(f)

        if len(self.sortmerna) == 0:
            log.debug("Could not find any SortMeRNA data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} logs".format(len(self.sortmerna)))
        self.write_data_file(self.sortmerna, 'multiqc_sortmerna')
        log.debug(self.sortmerna)
        # Add rRNA rate to the general stats table
        headers = OrderedDict()
        headers['rRNA_pct'] = {
            'title': '% rRNA',
            'description': '% rRNA',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'OrRd',
            'format': '{:.1f}%'
        }
        self.general_stats_addcols(self.sortmerna, headers)

        # Make barplot
        self.intro += self.sortmerna_overall_barplot()
        self.intro += self.sortmerna_detailed_barplot()

    def parse_sortmerna(self, f):
        s_name = None
        post_results_start = False
        post_database_start = False
        db_number = 0
        err = False

        for l in f["f"]:
            if "Reads file" in l:
                s_name = os.path.basename(l.split(" ")[-1]).strip().split(".")[0]
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

    def sortmerna_overall_barplot (self):
        keys = OrderedDict()
        keys["non_rRNA"] = { 'color': '#a6cee3', 'name': 'Other' }
        keys["rRNA"] = { 'color': '#e31a1c', 'name': 'rRNA' }
        pconfig = {
            'id': 'SortMeRNA overall',
            'title': 'SortMeRNA Other vs rRNA',
            'ylab': 'Reads'
        }
        log.debug(keys)
        return plots.bargraph.plot(self.sortmerna, keys, pconfig)

    def sortmerna_detailed_barplot (self):
        """ Make the HighCharts HTML to plot the sortmerna rates """

        colors = ["#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
                  "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
                  "#6a3d9a", "#ffff99", "#b15928"]

        # Specify the order of the different possible categories
        keys = OrderedDict()
        metrics = set()
        for sample in self.sortmerna:
            for key in self.sortmerna[sample]:
                if not key in ["total", "rRNA", "non_rRNA"] and not "_pct" in key:
                    metrics.add(key)

        col_index = 0
        for key in metrics:
            keys[key] = { 'color': colors[col_index], 'name': key.replace("_count","") }
            col_index = col_index + 1
        # Config for the plot
        pconfig = {
            'id': 'SortMeRNA detailed',
            'title': 'SortMeRNA hits',
            'ylab': 'Reads'
        }

        return plots.bargraph.plot(self.sortmerna, keys, pconfig)
