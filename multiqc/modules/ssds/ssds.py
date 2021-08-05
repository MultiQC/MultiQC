#!/usr/bin/env python

""" MultiQC module to parse output from SSDS reports """

from __future__ import print_function
from collections import OrderedDict

import logging
import re
import sys

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule

from multiqc.plots import linegraph, bargraph, scatter, table, heatmap, beeswarm

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """This MultiQC module supports output from two types of SSDS report.
    Parse SSDS report: Output from parsing an SSDS BAM file.
    SSDS SPoT report: Output of assessing SPoT of SSDS in genomic intervals"""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="SSDS QC reports",
            anchor="ssds",
            target="SSDS QC",
            href="http://genome.cshlp.org/content/early/2012/03/20/gr.130583.111.full.pdf",
            info=" Statistics for Single Stranded DNA Sequencing (SSDS).",
        )

        # Set up class objects to hold parsed data
        self.sections = list()
        self.general_stats_headers = OrderedDict()
        self.general_stats_data = dict()
        n = dict()

        # Call submodule functions
        n["ssds"] = self.parse_reports()
        if n["ssds"] > 0:
            log.info("Found {} parse SSDS reports".format(n["ssds"]))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)

    def parse_reports(self):
        """Find ssds stats logs and parse their data"""

        # Create dictionaries for stats and histogram data
        self.ssds_stats = dict()
        self.histograms = dict()
        self.heatMapList = dict()
        self.SPoT_headers = OrderedDict()

        self.heatXcats = []
        self.heatYcats = []

        # Loop through all files in folder for parse SSDS logs
        for f in self.find_log_files("ssds/parse"):

            # Chop off the file extension (save space on screen)
            reportName = f["s_name"]
            reportName = re.sub("\..+", "", reportName)

            parsed_data = dict()

            # Loop through the file by line
            for line in f["f"].splitlines():

                #  Split line into columns
                sections = line.split("\t")

                # Lines starting with totinfo are general information
                # This excludes those lines for histogram data (ITR/uH Len ...etc).
                if not line.startswith("totinfo"):

                    ## Build histograms
                    myType = sections[0]
                    nlen = sections[1]
                    nCount = sections[2]

                    # Create a dictionary for this type unless one already exists
                    if not (myType in self.histograms):
                        self.histograms[myType] = dict()

                    # Create a dictionary for this file and type unless one already exists
                    if not (f["fn"] in self.histograms[myType]):
                        self.histograms[myType][f["fn"]] = dict()

                    # Add this value : count pair to the histogram
                    self.histograms[myType][f["fn"]][int(nlen)] = int(nCount)

                    continue

                ## Get Simple Alignment data
                field = sections[1].strip()
                field = field.replace(" ", "_")
                value = float(sections[2].strip())
                parsed_data[field] = value

            # If we got data from a parse_ssds report
            if len(parsed_data) > 0:

                # Calculate some customized values for our data
                parsed_data["Total_Aligned"] = (
                    parsed_data["ssDNA_fragments"]
                    + parsed_data["ssDNA_type2_fragments"]
                    + parsed_data["dsDNA_hiconf_fragments"]
                    + parsed_data["dsDNA_loconf_fragments"]
                    + parsed_data["unclassified_fragments"]
                )

                parsed_data["other"] = parsed_data["total_fragments"] - parsed_data["valid_fragments"]

                if parsed_data["total_fragments"] > 0:
                    parsed_data["Aligned_percent"] = (
                        parsed_data["Total_Aligned"] / parsed_data["total_fragments"]
                    ) * 100
                else:
                    parsed_data["Aligned_percent"] = 0
                    parsed_data["other"] = 0

                # Work out some percentages for aligned reads
                if "total_fragments" in parsed_data:

                    # Scan through all keys
                    for k in list(parsed_data.keys()):

                        # Exclude things that have already been calculated
                        if (
                            k != "total_fragments"
                            and k != "Total_Aligned"
                            and k != "Aligned_percent"
                            and k != "ssDNA_NR"
                            and k != "other"
                            and k != "adapter"
                            and k != "adapter_percent"
                        ):
                            # Add a percentage for all T1/T2/ds/unc
                            if parsed_data["Total_Aligned"] > 0:
                                parsed_data["{}_percent".format(k)] = (
                                    (parsed_data[k] + 1) / parsed_data["Total_Aligned"]
                                ) * 100
                            else:
                                parsed_data["{}_percent".format(k)] = 0

                # Overwrite duplicate
                if f["s_name"] in self.ssds_stats:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))

                self.add_data_source(f)
                self.ssds_stats[reportName] = parsed_data

        for f in self.find_log_files("ssds/spot"):

            # Chop off the file extension (save space on screen)
            reportName = f["s_name"]
            reportName = re.sub("\..+", "", reportName)

            SPoT_data = dict()

            # Loop through the file by line
            for line in f["f"].splitlines():

                #  Split line into columns
                sections = line.split("\t")

                if sections[0].endswith("_SPoT"):
                    reptype = "SPoT"
                    read_dets = sections[0].replace("_SPoT", "")
                    interval_dets = sections[1]
                    SPoT = sections[2]

                    if not (read_dets in self.heatMapList):
                        self.heatMapList[read_dets] = dict()

                    if not (reportName in self.heatMapList[read_dets]):
                        self.heatMapList[read_dets][reportName] = dict()

                    if not (interval_dets in self.heatMapList[read_dets][reportName]):
                        self.heatMapList[read_dets][reportName][interval_dets] = dict()

                    self.SPoT_headers[interval_dets] = {
                        "title": interval_dets,
                        "description": "fragments in intervals / total : " + interval_dets,
                        "max": 1,
                        "min": 0,
                        "suffix": "",
                        "scale": "RdYlGn",
                        "format": "{:.3f}",
                    }

                    if float(SPoT) < 0.0005:
                        self.heatMapList[read_dets][reportName][interval_dets] = "<0.0001"
                    else:
                        self.heatMapList[read_dets][reportName][interval_dets] = SPoT

                    self.heatXcats.append(interval_dets)
                    self.heatYcats.append(reportName)

        # If we have some data to show, carry on
        if len(self.ssds_stats) > 0:

            # Write parsed report data to a file
            self.write_data_file(self.ssds_stats, "multiqc_ssds_stats")

            ######################################### General Stats Table
            # This is where we populate the general statistics table
            # Each object is a columnar entry
            self.general_stats_headers["Total_SSDS_fragments"] = {
                "title": "Tot_frags",
                "description": "Count of SSDS read pairs in raw aligned BAM (millions)",
                "min": 0,
                "max": 1000,
                "scale": "PuBu",
                "modify": lambda x: x / 1000000,
            }

            self.general_stats_headers["Valid_SSDS_fragments"] = {
                "title": "Valid_frags",
                "description": "Count of valid SSDS fragments in raw aligned BAM (millions)",
                "min": 0,
                "max": 1000,
                "scale": "PuBu",
                "modify": lambda x: x / 1000000,
            }

            self.general_stats_headers["Valid_percent"] = {
                "title": "Valid (%)",
                "description": "Valid Fragments / Total pairs(%)",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:.2f}%",
            }

            self.general_stats_headers["ssDNA_fragments_percent"] = {
                "title": "ssDNA(%)",
                "description": "ssDNA (%)",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:.2f}%",
            }

            self.general_stats_headers["ssDNA_type2_fragments_percent"] = {
                "title": "Type-2 (%)",
                "description": "Type 2 ssDNA (%)",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:.2f}%",
            }

            self.general_stats_headers["dsDNA_hiconf_fragments_percent"] = {
                "title": "dsDNA(hi)(%)",
                "description": "dsDNA hi-confidence (%)",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:.2f}%",
            }

            self.general_stats_headers["dsDNA_loconf_fragments_percent"] = {
                "title": "dsDNA(low)(%)",
                "description": "dsDNA low-confidence (%)",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:.2f}%",
            }

            self.general_stats_headers["unclassified_fragments_percent"] = {
                "title": "Unc (%)",
                "description": "Unclassified  (%)",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:.2f}%",
            }

            ## Don't know what this does, but it does no harm
            for s_name in self.ssds_stats:
                if s_name not in self.general_stats_data:
                    self.general_stats_data[s_name] = dict()
                self.general_stats_data[s_name].update(self.ssds_stats[s_name])

            ######################################### Begin making sub-plots

            reads = {"min": 0, "modify": lambda x: float(x) / 1000000.0, "decimalPlaces": 2, "shared_key": "read_count"}

            ##################################################################################
            ## PLOT 1: Total alignment barplot
            # Specify the order and colors of the different types of aligned / unaligned reads
            b2keys = OrderedDict()
            b2keys["ssDNA_fragments"] = {"color": "#437bb1", "name": "ssDNA"}
            b2keys["ssDNA_type2_fragments"] = {"color": "#f7a35c", "name": "ssDNA (type 2)"}
            b2keys["dsDNA_hiconf_fragments"] = {"color": "#11a400", "name": "dsDNA (higher confidence)"}
            b2keys["dsDNA_loconf_fragments"] = {"color": "#0b6c00", "name": "dsDNA (lower confidence)"}
            b2keys["unclassified_fragments"] = {"color": "#b1084c", "name": "Unclassified"}
            b2keys["adapter"] = {"color": "#7fffd4", "name": "Adapter"}
            b2keys["other"] = {"color": "#696969", "name": "Other"}

            # Configure the total alignment plot
            b2config = {
                "id": "Total_Alignment_Stats_Plot",
                "title": "SSDS: Read-pair_types",
                "ylab": "# Read-pairs",
                "cpswitch_counts_label": "Number of Read-pairs",
            }

            # Add the Total alignment stats barplot to the page
            self.add_section(
                plot=bargraph.plot(self.ssds_stats, b2keys, b2config),
                name="Read pair alignment stats",
                description='<p>This plot shows the fraction of read-pairs derived from different sources. It includes unmapped / low quality / supplementary read-pairs in the "other" category.</p>',
                content='<p>This module parses the results from <a href="https://github.com/kevbrick/ssds_pipeline_accessory_scripts.git"><code>parse_SSDS_BAM.py</code> </a></p>',
            )

            ##################################################################################
            ## PLOT 2: SSDS alignment barplot
            # Specify the order and colors of the different types of aligned reads
            keys = OrderedDict()
            keys["ssDNA_fragments"] = {"color": "#437bb1", "name": "ssDNA"}
            keys["ssDNA_type2_fragments"] = {"color": "#f7a35c", "name": "ssDNA (type 2)"}
            keys["dsDNA_hiconf_fragments"] = {"color": "#11a400", "name": "dsDNA (higher confidence)"}
            keys["dsDNA_loconf_fragments"] = {"color": "#0b6c00", "name": "dsDNA (lower confidence)"}
            keys["unclassified_fragments"] = {"color": "#b1084c", "name": "Unclassified"}

            # Configure the SSDS Alignment barplot
            pconfig = {
                "id": "SSDS_Alignment_Stats_Plot",
                "title": "SSDS: Fragment_types",
                "ylab": "# Fragments",
                "save_file": False,
                "cpswitch_counts_label": "Number of Fragments",
            }

            # Add the SSDS barplot to the page
            self.add_section(
                plot=bargraph.plot(self.ssds_stats, keys, pconfig),
                name="SSDS alignment stats (excluding unaligned)",
                anchor="ssds-alignment-barplot",
                description="<p>This plot shows the fraction of fragments derived from different sources.</p>",
                content='<p>This module parses the results from <a href="https://github.com/kevbrick/ssds_pipeline_accessory_scripts.git"><code>parse_SSDS_BAM.py</code> </a></p>',
            )

            ##################################################################################
            ## PLOT 3: SPoT Heatmap
            # Specify the order and colors of the different types of aligned / unaligned reads
            self.SPoT_headers_sorted = sorted(self.SPoT_headers)

            # Add the SPoT stats to the page
            for dna in ["ssDNA", "ssDNA_type2", "dsDNA_hiconf", "dsDNA_loconf", "unclassified"]:
                # Configure the SPoT table
                pconfig = {
                    "id": "SPoT_stats_table_" + dna,
                    "table_title": "SPoT Statistics for " + dna,
                    "save_file": False,
                    "raw_data_fn": "multiqc_SPoT_stats",
                }
                self.add_section(
                    plot=table.plot(self.heatMapList[dna], self.SPoT_headers, pconfig),
                    name="SPoTs (" + dna + ")",
                    anchor="SPoT-stats-" + dna,
                    description="<p>The signal portion of tags (SPoT) is the fraction of reads that are found at a given set of genomic intervals. \
                                    Values range from 0 (no reads) to 1 (all reads). The interval names have been defined by the user. Intervals \
                                    followed by (R) indicate the SPoT at a randomly shuffled set of the same intervals - this gives a ballpark idea \
                                    of the expected overlap if reads are not enriched at these loci. </p>",
                    content='<p>This module parses the results from <a href="https://github.com/kevbrick/ssds_pipeline_accessory_scripts.git"><code>calculate_SPoT.py</code> </a></p>',
                )

            ##################################################################################
            ## PLOTs 4-20: SSDS length histograms

            # Loop through all 20 combinations of parameter and type
            # (ss,t2,dH,dL,un) X (frag,ITR,uH,Fillin)
            # hard coding the list is the easiest option here (no need to be clever)
            plot_descriptions = {
                "Fragment": "NOTE: SSDS fragments tend to be rather short < 100bp.",
                "ITR": 'SSDS read pairs are characterized by Inverted Terminal Repeats (ITRs) (same sequence at the 5 prime end of read 1 and the 3 prime end of read 2). ITR length can be modulated by changing the temperature of the SSDS protocol (see <a href="https://genome.cshlp.org/content/22/5/957.long">Khil et al., Genome Res.2012</a>)',
                "uH": "ITR formation is thought to be modulated by short microhomologies (uH) in the genome.",
                "FillIn": "The fill-in region indicates the region at the 3 prime end of that is the result of sythesis using the 5 prime end of read 1 as a template. The fill-in at the end of read 2 is not found in the reference genome, but matches the 5 prime end of read 1.",
            }

            for dna_type in ["ssDNA", "ssDNA_type2", "dsDNA_hiconf", "dsDNA_loconf", "unclassified"]:
                for val_type in ["Fragment", "ITR", "uH", "FillIn"]:
                    combo_type = dna_type + "_" + val_type
                    # If histogram dictionary exists for this combo
                    if len(self.histograms[combo_type]) > 0:

                        ## Make a percentage normalised version of the data
                        data_percent = {}

                        # Loop through key, value pairs for this histogram
                        for s_name, data in self.histograms[combo_type].items():

                            # initialize total to 0 (not sure why I do that)
                            total = 0

                            # Create percentage dictionary
                            data_percent[s_name] = OrderedDict()

                            # Get total for this histogram
                            total = float(sum(data.values()))

                            # Calculate percentages for this histogram
                            for k, v in data.items():
                                if v > 0:
                                    data_percent[s_name][k] = (v / total) * 100
                                else:
                                    data_percent[s_name][k] = 0

                        fig_label = combo_type.replace("_", " ")
                        # Configure histogram plot
                        histConfig = {
                            "id": combo_type + "_length",
                            "title": "SSDS: " + fig_label + " length",
                            "ylab": "Fragments (%)",
                            "save_file": False,
                            "xlab": fig_label + " length (bp)",
                            "xDecimals": False,
                            "ymin": 0,
                            "data_labels": [
                                {"name": "Percent", "ylab": "Fragments (%)"},
                                {"name": "Counts", "ylab": "Fragments (count)"},
                            ],
                        }

                        # Add histogram to multi-QC page
                        self.add_section(
                            plot=linegraph.plot([data_percent, self.histograms[combo_type]], histConfig),
                            name=fig_label,
                            anchor="ssds-stats" + combo_type,
                            description="<p>This plot shows the length distribution of SSDS "
                            + val_type
                            + "s for "
                            + dna_type
                            + " fragments. "
                            + plot_descriptions[val_type]
                            + ".</p>",
                            content='<p>This module parses the results from <a href="https://github.com/kevbrick/ssds_pipeline_accessory_scripts.git"><code>parse_SSDS_BAM.py</code> </a></p>',
                        )

        # Return the number of logs that were found
        return len(self.ssds_stats)
