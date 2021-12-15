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
    SSDS details report: Output from parsing an SSDS BAM file.
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
        n["SSDS_dets"] = self.parse_details_reports()
        if n["SSDS_dets"] > 0:
            log.info("Found {} SSDS details reports".format(n["SSDS_dets"]))

        # Call submodule functions
        n["SSDS_SPoT"] = self.parse_spot_reports()
        if n["SSDS_SPoT"] > 0:
            log.info("Found {} SSDS SPoT reports".format(n["SSDS_SPoT"]))

        # Exit if we didn't find anything
        if sum(n.values()) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        # Add to the General Stats table (has to be called once per MultiQC module)
        self.general_stats_addcols(self.general_stats_data, self.general_stats_headers)

        if n["SSDS_dets"] > 0:
            self.tot_alignment_plot()
            self.ssds_alignment_plot()
            self.itr_properties_plots()

        if n["SSDS_SPoT"] > 0:
            self.ssds_spot_heatmap()

    def parse_details_reports(self):

        ## Find ssds-parse logs and get data

        # Create dictionaries for stats and histogram data
        self.ssds_stats = dict()
        self.histograms = dict()

        # Loop through all files in folder for parse SSDS logs
        for f in self.find_log_files("ssds/parse"):

            ## Prevent loading of empty reports
            reportOK = False

            # Chop off the file extension (save space on screen)
            reportName = f["s_name"]
            reportName = re.sub("\..+", "", reportName)

            parsed_data = dict()

            # Loop through the file by line
            for line in f["f"].splitlines():
                reportOK = True

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
                if field == "filtered_fragments":
                    parsed_data["other"] = parsed_data[field]

            if reportOK:
                self.add_data_source(f)
                self.ssds_stats[reportName] = parsed_data

        # If we have some data to show, carry on
        if len(self.ssds_stats) > 0:

            # Write parsed report data to a file
            self.write_data_file(self.ssds_stats, "multiqc_ssds_stats")

            ######################################### General Stats Table
            # This is where we populate the general statistics table
            # Each object is a columnar entry
            self.general_stats_headers["ssDNA_fragments"] = {
                "title": "ssDNA",
                "description": "Count of SSDS read pairs in raw aligned BAM (millions)",
                "scale": "PuBu",
                "modify": lambda x: x / 1000000,
            }

            for s_name in self.ssds_stats:
                if s_name not in self.general_stats_data:
                    self.general_stats_data[s_name] = dict()
                self.general_stats_data[s_name].update(self.ssds_stats[s_name])

        # Return the number of logs that were found
        return len(self.ssds_stats)

    def parse_spot_reports(self):

        ## Find SSDS SPot report logs and parse their data

        # Create dictionaries
        self.SPoT_values = dict()
        self.ssds_SPoT_stats = dict()

        for f in self.find_log_files("ssds/spot"):

            # Chop off the file extension (save space on screen)
            reportName = f["s_name"]
            reportName = re.sub("\..+", "", reportName)

            reportOK = False

            # Loop through the file by line
            for line in f["f"].splitlines():

                #  Split line into columns
                sections = line.split("\t")

                if sections[0].endswith("_SPoT"):
                    reportOK = True
                    reptype = "SPoT"
                    read_dets = sections[0].replace("_SPoT", "")
                    interval_dets = sections[1]
                    SPoT = sections[2]

                    if not (read_dets in self.SPoT_values):
                        self.SPoT_values[read_dets] = dict()

                    if not (reportName in self.SPoT_values[read_dets]):
                        self.SPoT_values[read_dets][reportName] = dict()

                    if not (interval_dets in self.SPoT_values[read_dets][reportName]):
                        self.SPoT_values[read_dets][reportName][interval_dets] = dict()

                    if float(SPoT) < 0.0005:
                        self.SPoT_values[read_dets][reportName][interval_dets] = "0"
                    else:
                        self.SPoT_values[read_dets][reportName][interval_dets] = float(SPoT) * 100

            if reportOK:
                self.add_data_source(f)
                self.ssds_SPoT_stats[reportName] = self.SPoT_values[read_dets]

        # Return the number of logs that were found
        return len(self.ssds_SPoT_stats)

    def tot_alignment_plot(self):

        ## PLOT 1: Total alignment barplot

        sample_color_keys = OrderedDict()
        sample_color_keys["ssDNA_fragments"] = {"color": "#437bb1", "name": "ssDNA"}
        sample_color_keys["ssLow_fragments"] = {"color": "#103150", "name": "ssDNA (low confidence)"}
        sample_color_keys["type2_fragments"] = {"color": "#f7a35c", "name": "Ambiguous (ss/ds DNA)"}
        sample_color_keys["dsDNA_fragments"] = {"color": "#11a400", "name": "dsDNA"}
        sample_color_keys["unclassified_fragments"] = {"color": "#b1084c", "name": "Unclassified"}
        sample_color_keys["adapter"] = {"color": "#7fffd4", "name": "Adapter"}
        sample_color_keys["other"] = {"color": "#696969", "name": "Other"}

        # Configure the total alignment plot
        pconfig = {
            "id": "Total_Alignment_Stats_Plot",
            "title": "SSDS: Read-pair_types",
            "ylab": "# Read-pairs",
            "cpswitch_counts_label": "Number of Read-pairs",
        }

        # Add the Total alignment stats barplot to the page
        self.add_section(
            plot=bargraph.plot(self.ssds_stats, sample_color_keys, pconfig),
            name="Read pair alignment stats",
            description='<p>This plot shows the fraction of read-pairs derived from different sources. It includes unmapped / low quality / supplementary read-pairs in the "other" category.</p>',
            content='<p>This module parses the results from <a href="https://github.com/kevbrick/ssds_pipeline_accessory_scripts.git"><code>parse_SSDS_BAM.py</code> </a></p>',
        )

    def ssds_alignment_plot(self):

        ## PLOT 2: SSDS alignment barplot

        keys = OrderedDict()
        keys["ssDNA_fragments"] = {"color": "#437bb1", "name": "ssDNA"}
        keys["ssLow_fragments"] = {"color": "#103150", "name": "ssDNA (low confidence)"}
        keys["type2_fragments"] = {"color": "#f7a35c", "name": "Ambiguous (ss/ds DNA)"}
        keys["dsDNA_fragments"] = {"color": "#11a400", "name": "dsDNA"}
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

    def itr_properties_plots(self):

        ## Plot 4: SSDS properties / Length histograms
        plot_description = """<table style="width:100%">
            <tr style="vertical-align: top; text-align: left"><td>Fragment length:</td><td>Also known as "insert size". Fragments from SSDS experiments are generally short (<100bp).</td>
            <tr style="vertical-align: top; text-align: left"><td>Inverted Terminal Repeats (ITRs):</td><td>SSDS read pairs are characterized by Inverted Terminal Repeats (ITRs) (the same sequence at the 5 prime end of read 1 and the 3 prime end of read 2).<br>ITRs with a fill-in are characteristic of ssDNA (see <a href="https://genome.cshlp.org/content/22/5/957.long">Khil et al., Genome Res.2012</a>).</td>
            <tr style="vertical-align: top; text-align: left"><td>ITR micro-homology (uH):</td><td>ITR formation is modulated by short microhomologies (uH) in the genome.</td>
            <tr style="vertical-align: top; text-align: left"><td>ITR Fill-In:</td><td>The fill-in region arises during end-repair; DNA synthesis extends the 3' end of ssDNA using the 5 prime end of the fragment as a template.<br>The fill-in occurs at the 3 prime end of read 2 is not found in the reference genome, but matches the 5 prime end of read 2 and is not found in the genome.</td></tr></table>"""

        plot_descriptions = {}
        plot_descriptions["ssDNA"] = (
            "Fragments designated <b>ssDNA</b> are derived from single-stranded DNA.<hr>"
            + plot_description
            + "<hr><i>Selection criteria:</i> <b>ITR > 5 bp</b> & <b>Fill-in > 2 bp</b> & <b>uHomology > 0 bp</b><hr>"
        )
        plot_descriptions["ssLow"] = (
            "Fragments designated as <b>ssLow</b> are likely derived from single-stranded DNA. However, this is a low-confidence designation as these fragments may also be derived from dsDNA. Not routinely used in ssDNA-based analyses.<hr>"
            + plot_description
            + "<hr><i>Selection criteria:</i> <b>Fill-in > 0 bp</b> but <b>not classified as ssDNA</b><hr>"
        )
        plot_descriptions["type2"] = (
            "Fragments designated as <b>type2</b> may be derived from ssDNA or dsDNA. These fragments are not routinely used in ssDNA-based analyses..<hr>"
            + plot_description
            + "<hr><i>Selection criteria:</i> <b>ITR > 3 bp</b> & <b>NO fill-in</b><hr>"
        )
        plot_descriptions["dsDNA"] = (
            "Fragments designated as <b>dsDNA</b> are likely derived from double-stranded DNA. It remains possible that these fragments are derived from single-stranded DNA, particularly in libraries with ssDNA enrichment.<hr>"
            + plot_description
            + "<hr><i>Selection criteria:</i> <b>ITR < 3 bp</b> & <b>NO fill-in</b><hr>"
        )
        plot_descriptions["unclassified"] = (
            "Fragments that cannot be classified based on the ITR structure. Designated as <b>unclassified</b>.<hr>"
            + plot_description
            + "<hr><i>Selection criteria:</i> <b>All remaining reads</b> not classified by other criteria.<hr>"
        )

        for dna_type in ["ssDNA", "ssLow", "type2", "dsDNA", "unclassified"]:

            ok_to_plot = False
            ## Make a percentage normalised version of the data
            data_count = {}
            data_percent = {}

            for prop_type in ["ITR", "uH", "FillIn", "Fragment"]:
                combo_type = dna_type + "_" + prop_type

                data_count[prop_type] = {}
                data_percent[prop_type] = {}

                # If histogram dictionary exists for this combo
                if combo_type in self.histograms.keys():
                    if len(self.histograms[combo_type]) > 0:
                        ok_to_plot = True

                        # Loop through key, value pairs for this histogram
                        for s_name, data in self.histograms[combo_type].items():
                            maxXval = {}
                            maxXval[prop_type] = -1

                            total = 0
                            # Create percentage dictionary
                            data_count[prop_type][s_name] = {}
                            data_percent[prop_type][s_name] = {}

                            # Get total for this histogram
                            total = float(sum(data.values()))

                            # Calculate percentages for this histogram
                            for k, v in data.items():
                                if (maxXval[prop_type] + 1) < k:
                                    for l in range(maxXval[prop_type] + 1, k):
                                        data_count[prop_type][s_name][l] = 0
                                        data_percent[prop_type][s_name][l] = 0

                                data_count[prop_type][s_name][k] = v
                                if v > 0:
                                    data_percent[prop_type][s_name][k] = (v / total) * 100
                                else:
                                    data_percent[prop_type][s_name][k] = 0
                                if k > maxXval[prop_type]:
                                    maxXval[prop_type] = k

            if ok_to_plot:
                # Configure histogram plot
                histConfig = {
                    "id": dna_type + "_fragment_properties",
                    "title": "Module SSDS: " + dna_type + " fragment properties",
                    "save_file": True,
                    "xDecimals": False,
                    "ymin": 0,
                    "ylab": "Fragments",
                    "categories": True,
                    "data_labels": [
                        {"name": "ITR (Count)", "ylab": "Fragments (#)", "xlab": "Total ITR length (nt)"},
                        {"name": "ITR (%)", "ylab": "Fragments (%)", "xlab": "Total ITR length (nt)"},
                        {
                            "name": "Micro-homology (Count)",
                            "ylab": "Fragments (#)",
                            "xlab": "ITR microhomology length (nt)",
                        },
                        {
                            "name": "Micro-homology (%)",
                            "ylab": "Fragments (%)",
                            "xlab": "ITR microhomology length (nt)",
                        },
                        {"name": "Fill-in (Count)", "ylab": "Fragments (#)", "xlab": "ITR fill-in length (nt)"},
                        {"name": "Fill-in (%)", "ylab": "Fragments (%)", "xlab": "ITR fill-in length (nt)"},
                        {"name": "Fragment (Count)", "ylab": "Fragments (#)", "xlab": "Fragment length (nt)"},
                        {"name": "Fragment (%)", "ylab": "Fragments (%)", "xlab": "Fragment length (nt)"},
                    ],
                }

                plot_data = [
                    data_count["ITR"],
                    data_percent["ITR"],
                    data_count["uH"],
                    data_percent["uH"],
                    data_count["FillIn"],
                    data_percent["FillIn"],
                    data_count["Fragment"],
                    data_percent["Fragment"],
                ]

                # KB: 10/08/21: Adding this linegraph is slow.
                # Add histogram to multi-QC page
                self.add_section(
                    plot=linegraph.plot(plot_data, histConfig),
                    name="Frag. props: " + dna_type,
                    anchor="ssds-stats" + combo_type,
                    description="<p>This plot shows the length distributions for fragment properties for "
                    + dna_type
                    + " fragments from SSDS.<br><br>"
                    + plot_descriptions[dna_type]
                    + "</p>",
                    content='<p>This module parses the results from <a href="https://github.com/kevbrick/ssds_pipeline_accessory_scripts.git"><code>parse_SSDS_BAM.py</code> </a></p>',
                )

    def ssds_spot_heatmap(self):

        ## PLOT 3: Heatmap showing SPoT breakdown by type for every sample

        data = []
        spot_vals = OrderedDict()

        help_text = {}
        help_text[
            "ssDNA"
        ] = """<tr style="vertical-align: top; text-align: left"><td><b>ss</b> (ssDNA)        : Fragments unambiguously derived from single-stranded DNA (ssDNA)</td>"""
        help_text[
            "ssLow"
        ] = """<tr style="vertical-align: top; text-align: left"><td><b>sl</b> (ssLow)        : Low confidence ssDNA. With a low likelihood, these fragments may also be derived from dsDNA. Not routinely used in ssDNA-based analyses.</td>"""
        help_text[
            "type2"
        ] = """<tr style="vertical-align: top; text-align: left"><td><b>t2</b> (type2)        : Fragments may be derived from ss- or ds-DNA. Not routinely used in ssDNA-based analyses.</td>"""
        help_text[
            "dsDNA"
        ] = """<tr style="vertical-align: top; text-align: left"><td><b>ds</b> (dsDNA)        : Fragments likely derived from double-stranded DNA. It remains possible that these fragments are derived from single-stranded DNA, particularly in libraries with ssDNA enrichment. </td>"""
        help_text[
            "unclassified"
        ] = """<tr style="vertical-align: top; text-align: left"><td><b>un</b> (unclassified) : Ambiguous, unclassified fragments.</td>"""

        hm = self.SPoT_values

        # dna_types = self.SPoT_values.keys()
        dna_types = []
        help_dets = ""
        for d in ["ssDNA", "ssLow", "type2", "dsDNA", "unclassified"]:
            if d in self.SPoT_values.keys():
                dna_types.append(d)
                help_dets = help_dets + help_text[d]

        short_dna_type = OrderedDict()
        short_dna_type["ssDNA"] = "ss"
        short_dna_type["ssLow"] = "sl"
        short_dna_type["type2"] = "t2"
        short_dna_type["dsDNA"] = "ds"
        short_dna_type["unclassified"] = "un"

        sample_names = sorted(self.SPoT_values["ssDNA"].keys())

        interval_names = []
        for s in sample_names:
            for k in sorted(self.SPoT_values["ssDNA"][s]):
                if k not in interval_names:
                    interval_names.append(k)

        s_names = []

        for d in dna_types:
            for s in sample_names:
                s_names.append("(" + short_dna_type[d] + ")" + s)

        for d in dna_types:
            for s in sample_names:
                row = []
                for i in interval_names:
                    try:
                        row.append(float(self.SPoT_values[d][s][i]))
                    except KeyError:
                        row.append(0)
                data.append(row)

        pconfig = {
            "id": "ssds-spot-heatmap",
            "title": "SSDS: Signal Percentage of Tags (%)",
            "xTitle": "Sample",
            "yTitle": "Interval",
            "square": False,
            "colstops": [
                [0, "#ffffff"],
                [0.001, "#fefce9"],
                [0.50, "#ffc265"],
                [1.00, "#ff6262"],
            ],
            "decimalPlaces": 0,
            "legend": False,
            "datalabels": True,
            "xcats_samples": False,
            "borderWidth": 1,
        }

        self.add_section(
            name="SSDS SPoTs",
            anchor="ssds_spot_heatmap",
            description="""
                Signal Percentage of Tags (SPoT) for all samples (%). Colors indicate the value (0 / no data =white; 
                Otherwise, increasing SPoT from yellow to orange to red). Intervals annotated as (R) represent
                the SPoT when the intervals are randomly shuffled in the genome (bedtools shuffle -chrom). This 
                provides a naive, but useful estimate of random expectation for a non-enriched library.<br><br>
                <table style="width:100%">"""
            + help_dets
            + "</table><br>",
            helptext="""
                The Signal Percentage of Tags (SPoT) represents the percentage of sequencing reads found in
                a set of genomic intervals. Higher numbers indicate that the library was enriched for reads
                in that location. The SSDS report may also contain intervals annotated as (R); these represent
                the SPoT when the intervals are randomly shuffled in the genome (bedtools shuffle -chrom). This
                represents a reasonable expectation of random overlap, however this very simple estimate should
                be formally validated more robustly.  
            """,
            plot=heatmap.plot(data, interval_names, s_names, pconfig),
        )
