""" MultiQC module to parse output from SortMeRNA """


import logging
import os
import re
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="SortMeRNA",
            anchor="sortmerna",
            href="http://bioinfo.lifl.fr/RNA/sortmerna/",
            info="is a program for filtering, mapping and OTU-picking NGS reads in metatranscriptomic and metagenomic data.",
            doi="10.1093/bioinformatics/bts611",
        )

        # Parse logs
        self.sortmerna = dict()
        for f in self.find_log_files("sortmerna", filehandles=True):
            self.parse_sortmerna(f)
            self.add_data_source(f)

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, f["s_name"])

        # Filter to strip out ignored sample names
        self.sortmerna = self.ignore_samples(self.sortmerna)

        if len(self.sortmerna) == 0:
            raise ModuleNoSamplesFound

        log.info("Found {} logs".format(len(self.sortmerna)))
        self.write_data_file(self.sortmerna, "multiqc_sortmerna")

        # Get custom table header, default to 'rRNA'
        tab_header = getattr(config, "sortmerna", {}).get("seqname", "% rRNA")

        # Add rRNA rate to the general stats table
        headers = OrderedDict()
        headers["rRNA_pct"] = {
            "title": tab_header,
            "description": "Percentage of reads matched to a SortMeRNA database",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "OrRd",
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
                parts = re.split(":|=", l)
                s_name = self.clean_s_name(parts[-1], f)
                self.sortmerna[s_name] = dict()
            if "Results:" in l and not post_results_start:  # old versions
                post_results_start = True
            if "Total reads = " in l and not post_results_start:  # v4.2.0 onwards
                post_results_start = True
            if not post_results_start:
                continue
            if post_results_start and not post_database_start:
                if "Total reads =" in l:
                    m = re.search("\d+", l)
                    if m:
                        self.sortmerna[s_name]["total"] = int(m.group())
                    else:
                        err = True
                elif "Total reads passing" in l:
                    m = re.search("\d+", l)
                    if m:
                        self.sortmerna[s_name]["rRNA"] = int(m.group())
                        self.sortmerna[s_name]["rRNA_pct"] = (
                            float(self.sortmerna[s_name]["rRNA"]) / float(self.sortmerna[s_name]["total"]) * 100
                        )
                    else:
                        err = True
                elif "Total reads failing" in l:
                    m = re.search("\d+", l)
                    if m:
                        self.sortmerna[s_name]["non_rRNA"] = int(m.group())
                        self.sortmerna[s_name]["non_rRNA_pct"] = (
                            float(self.sortmerna[s_name]["non_rRNA"]) / float(self.sortmerna[s_name]["total"]) * 100
                        )
                    else:
                        err = True
            if post_database_start:
                # Stop when we hit empty lines
                if not l.strip():
                    break

                db_number = db_number + 1
                parts = re.split("\t+", l.strip())
                if len(parts) == 2:
                    db = os.path.splitext(os.path.basename(parts[0]))[0]
                    pct = float(parts[1].replace("%", ""))
                    count = int(self.sortmerna[s_name]["total"]) * (pct / 100.0)
                    self.sortmerna[s_name][db + "_pct"] = pct
                    self.sortmerna[s_name][db + "_count"] = count
                else:
                    err = True

            if "By database:" in l or "Coverage by database:" in l:
                post_database_start = True
        if err:
            log.warning("Error parsing data in: " + s_name)
            self.sortmerna.pop(s_name, "None")
        s_name = None

    def sortmerna_detailed_barplot(self):
        """Make the HighCharts HTML to plot the sortmerna rates"""

        # Specify the order of the different possible categories
        keys = {}
        metrics = set()
        for sample in self.sortmerna:
            for key in self.sortmerna[sample]:
                if not key in ["total", "rRNA", "non_rRNA"] and not "_pct" in key:
                    metrics.add(key)

        for key in metrics:
            keys[key] = {"name": key.replace("_count", "")}

        # Config for the plot
        pconfig = {"id": "sortmerna-detailed-plot", "title": "SortMeRNA: Hit Counts", "ylab": "Reads"}

        self.add_section(plot=bargraph.plot(self.sortmerna, sorted(keys), pconfig))
