import json
from collections import OrderedDict

import numpy as np

from multiqc.plots import bargraph, linegraph, scatter, table
from multiqc.utils import report

"""
Functions for plotting per sample information of bases2fastq
"""


def tabulate_sample_stats(sample_data, group_lookup_dict, sample_color):
    """
    Tabulate general information and statistics per sample
    """
    plot_content = dict()
    for s_name in sample_data.keys():
        general_stats = dict()
        general_stats.update({"num_polonies_sample": sample_data[s_name]["NumPolonies"]})
        general_stats.update({"mean_base_quality": sample_data[s_name]["QualityScoreMean"]})
        general_stats.update({"yield_sample": sample_data[s_name]["Yield"]})
        general_stats.update({"percent_q30_sample": sample_data[s_name]["PercentQ30"]})
        general_stats.update({"percent_q40_sample": sample_data[s_name]["PercentQ40"]})
        general_stats.update({"group": group_lookup_dict[s_name]})
        plot_content.update({s_name: general_stats})

    headers = OrderedDict()

    headers["group"] = {
        "title": "Group",
        "description": "Assigned group name. The group name will be project name if given, or run name",
        "bgcols": sample_color,
    }
    headers["num_polonies_sample"] = {
        "title": "#Polonies",
        "description": "Number of polonies assigned to this sample",
        "format": "{d}",
        "min": 0,
        "scale": "Blues",
    }
    headers["yield_sample"] = {
        "title": "Yield(Gb)",
        "description": "Percent of reads with perfect index (0 mismatches)",
        "scale": "Greens",
    }
    headers["mean_base_quality"] = {
        "title": "Average Base Quality",
        "description": "Average base quality across R1/R2",
        "min": 0,
        "scale": "Spectral",
    }
    headers["percent_q30_sample"] = {
        "title": "% Bases Q30",
        "description": "Percent of reads with perfect index (0 mismatches)",
        "max": 100,
        "min": 0,
        "scale": "RdYlGn",
        "suffix": "%",
    }
    headers["percent_q40_sample"] = {
        "title": "% Bases Q40",
        "description": "Percent of reads with perfect index (0 mismatches)",
        "max": 100,
        "min": 0,
        "scale": "RdYlGn",
        "suffix": "%",
    }

    config = {"description": "Table of per sample key informations", "no_beeswarm": True}

    plot_name = "Sample QC metrics table"
    plot_html = table.plot(plot_content, headers, pconfig=config)
    anchor = "sample_qc_metrics_table"
    description = "table of general QC metrics by sample"
    helptext = """
    This section shows numbers of some metrics that indicate the quality of each bases. 
    
    Polony numbers: total number of polonies assigned to this sample. Each polony yields one read1 and one read2
    
    Yield: total volume of data that has been assigned to this sample

    Percentage Q30: percentage of bases that has >=30 base quality (error rate <= 10^-3) 
    
    Percentage Q40: percentage of bases that has >=40 base quality (error rate <= 10^-4)
    
    Average Base Quality: average base qualities across all bases
    
    Group: A group tag for assigning colors in the plot. The default group of each sample will be their sequencing
    run name. To customize group tags, you can either
    
    1) Set the project name when running bases2fastq. In this case the group tags will be project name.
    
    2) Generate a csv file that has the columns "Run Name","Sample Name" and "Group". An easy way of 
    generating this csv will be run multiqc with default grouping, copy this table into a csv file, 
    and then change group name manually. Then move the csv file into any location under the analyzed 
    directory.
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content


def sequence_content_plot(sample_data, group_lookup_dict, color_dict):
    """Create the epic HTML for the FastQC sequence content heatmap"""

    # Prep the data
    data = dict()

    r1r2_split = 0
    for s_name in sorted(sample_data.keys()):
        data[s_name] = {}
        R1 = sample_data[s_name]["Reads"][0]["Cycles"]
        r1r2_split = max(r1r2_split, len(R1))

    for s_name in sorted(sample_data.keys()):
        data[s_name] = {}
        R1 = sample_data[s_name]["Reads"][0]["Cycles"]
        for cycle in range(len(R1)):
            base_no = str(cycle + 1)
            data[s_name].update({base_no: {base.lower(): R1[cycle]["BaseComposition"][base] for base in "ACTG"}})
            data[s_name][base_no]["base"] = base_no
            tot = sum([data[s_name][base_no][base] for base in ["a", "c", "t", "g"]])
            for base in ["a", "c", "t", "g"]:
                data[s_name][base_no][base] = (float(data[s_name][base_no][base]) / float(tot)) * 100.0

        R2 = sample_data[s_name]["Reads"][1]["Cycles"]
        for cycle in range(len(R2)):
            base_no = str(cycle + 1 + r1r2_split)
            data[s_name].update({base_no: {base.lower(): R2[cycle]["BaseComposition"][base] for base in "ACTG"}})
            data[s_name][base_no]["base"] = base_no
            tot = sum([data[s_name][base_no][base] for base in ["a", "c", "t", "g"]])
            for base in ["a", "c", "t", "g"]:
                data[s_name][base_no][base] = (float(data[s_name][base_no][base]) / float(tot)) * 100.0
    html = """<div id="fastqc_per_base_sequence_content_plot_div">
        <div class="alert alert-info">
            <span class="glyphicon glyphicon-hand-up"></span>
            Click a sample row to see a line plot for that dataset.
        </div>
        <h5><span class="s_name text-primary"><span class="glyphicon glyphicon-info-sign"></span> Rollover for sample name</span></h5>
        <button id="fastqc_per_base_sequence_content_export_btn"><span class="glyphicon glyphicon-download-alt"></span> Export Plot</button>
        <div class="fastqc_seq_heatmap_key">
            Position: <span id="fastqc_seq_heatmap_key_pos">-</span>
            <div><span id="fastqc_seq_heatmap_key_t"> %T: <span>-</span></span></div>
            <div><span id="fastqc_seq_heatmap_key_c"> %C: <span>-</span></span></div>
            <div><span id="fastqc_seq_heatmap_key_a"> %A: <span>-</span></span></div>
            <div><span id="fastqc_seq_heatmap_key_g"> %G: <span>-</span></span></div>
        </div>
        <div id="fastqc_seq_heatmap_div" class="fastqc-overlay-plot">
            <div id="{id}" class="fastqc_per_base_sequence_content_plot hc-plot has-custom-export">
                <canvas id="fastqc_seq_heatmap" height="100%" width="800px" style="width:100%;"></canvas>
            </div>
        </div>
        <div class="clearfix"></div>
    </div>
    <script type="application/json" class="fastqc_seq_content">{d}</script>
    <script type="application/json" class="bases2fastq_group_color">{c}</script>
    <script> R1cyclenum = {r}</script>
    """.format(
        # Generate unique plot ID, needed in mqc_export_selectplots
        id=report.save_htmlid("per_base_composition_plot"),
        d=json.dumps(["bases2fastq", data]),
        c=json.dumps(["bases2fastq", color_dict]),
        r=r1r2_split,
    )
    plot_name = "Per Cycle Sequence Composition"
    anchor = "per_cycle_sequence_content"
    description = "The proportion of each base position for which each of the four normal DNA bases has been called."
    helptext = """
    The heatmap is formatted in a similar way as fastqc per base sequence content. The differences are:
    1) The pass/fail flag in the original fastqc module are replaced with a color reflecting the group tag of the samples
    2) Read1 and Read2 is co-plotted and separate by a white band in the heatmap and a red dashed in each individual sample plot.

    The following is a copy of helptext of the fastqc module:
    To enable multiple samples to be shown in a single plot, the base composition data
    is shown as a heatmap. The colours represent the balance between the four bases:
    an even distribution should give an even muddy brown colour. Hover over the plot
    to see the percentage of the four bases under the cursor.

    **To see the data as a line plot, as in the original FastQC graph, click on a sample track.**

    From the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html):

    _Per Base Sequence Content plots out the proportion of each base position in a
    file for which each of the four normal DNA bases has been called._

    _In a random library you would expect that there would be little to no difference
    between the different bases of a sequence run, so the lines in this plot should
    run parallel with each other. The relative amount of each base should reflect
    the overall amount of these bases in your genome, but in any case they should
    not be hugely imbalanced from each other._

    _It's worth noting that some types of library will always produce biased sequence
    composition, normally at the start of the read. Libraries produced by priming
    using random hexamers (including nearly all RNA-Seq libraries) and those which
    were fragmented using transposases inherit an intrinsic bias in the positions
    at which reads start. This bias does not concern an absolute sequence, but instead
    provides enrichement of a number of different K-mers at the 5' end of the reads.
    Whilst this is a true technical bias, it isn't something which can be corrected
    by trimming and in most cases doesn't seem to adversely affect the downstream
    analysis._
    """
    plot_content = data
    return html, plot_name, anchor, description, helptext, plot_content


def plot_per_cycle_N_content(sample_data, group_lookup_dict, color_dict):
    data = dict()
    r1r2_split = 0
    for s_name in sorted(sample_data.keys()):
        data[s_name] = {}
        R1 = sample_data[s_name]["Reads"][0]["Cycles"]
        R1_cycle_num = len(R1)
        r1r2_split = max(r1r2_split, R1_cycle_num)

    for s_name in sorted(sample_data.keys()):
        data[s_name] = {}
        R1 = sample_data[s_name]["Reads"][0]["Cycles"]
        R1_cycle_num = len(R1)
        for cycle in range(len(R1)):
            base_no = str(cycle + 1)
            if sum(R1[cycle]["BaseComposition"].values()) == 0:
                data[s_name].update({base_no: 0})
            else:
                data[s_name].update(
                    {base_no: R1[cycle]["BaseComposition"]["N"] / sum(R1[cycle]["BaseComposition"].values()) * 100.0}
                )

        R2 = sample_data[s_name]["Reads"][1]["Cycles"]
        R2_cycle_num = len(R2)
        for cycle in range(len(R2)):
            base_no = str(cycle + 1 + r1r2_split)
            if sum(R2[cycle]["BaseComposition"].values()) == 0:
                data[s_name].update({base_no: 0})
            else:
                data[s_name].update(
                    {base_no: R2[cycle]["BaseComposition"]["N"] / sum(R2[cycle]["BaseComposition"].values()) * 100.0}
                )

    plot_content = data
    config = {
        "xlab": "cycle",
        "ylab": "Base Quality",
        "ymax": 100,
        "xPlotLines": [{"color": "#FF0000", "width": 2, "value": r1r2_split, "dashStyle": "Dash"}],
        "colors": color_dict,
        "ymin": 0,
        "id": "per_cycle_n_content",
        "title": "bases2fastq: per cycle N content percentage",
        "ylab": "Percentage",
    }
    plot_html = linegraph.plot(plot_content, pconfig=config)
    plot_name = "Per Cycle N Content"
    anchor = "n_content"
    description = "N content by cycle"
    helptext = """
    This section plots the percentage of unidentified bases ("N" bases) by each sequencing cycle. 
    Read 1 and Read 2 are separated by a red dashed line
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content


def plot_per_read_gc_hist(sample_data, group_lookup_dict, sample_color):
    """
    Plot gc histogram per sample
    """
    gc_hist_dict = dict()
    for s_name in sample_data.keys():
        R1_gc_counts = sample_data[s_name]["Reads"][0]["PerReadGCCountHistogram"]
        R2_gc_counts = sample_data[s_name]["Reads"][1]["PerReadGCCountHistogram"]
        R1R2_gc_counts = [r1 + r2 for r1, r2 in zip(R1_gc_counts, R2_gc_counts)]
        totalReads = sum(R1R2_gc_counts)
        gc_hist_dict.update({s_name: {}})
        RLen = len(R1_gc_counts)
        for gc in range(0, RLen):
            gc_hist_dict[s_name].update({gc / RLen * 100: R1R2_gc_counts[gc] / totalReads * 100})

    # perReadQualityHistogram
    plot_content = gc_hist_dict
    plot_function = linegraph.plot

    config = {
        "description": "GC",
        "xlab": "GC content",
        "ylab": "Percentage",
        "colors": sample_color,
        "id": "gc_hist",
        "title": "bases2fastq: per sample GC content histogram",
        "ylab": "Percentage",
    }
    plot_name = "per sample GC histogram"
    plot_html = linegraph.plot(plot_content, pconfig=config)
    anchor = "gc_histogram"
    description = "Histogram of distributions of percentage GC in each read"
    helptext = """
    This section plots the distribution of percentage GC in each reads (range: 0-100)
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content


def plot_adapter_content(sample_data, group_lookup_dict, sample_color):
    """
    Plot adapter content per sample
    """
    plot_content = dict()

    r1r2_split = 0
    for s_name in sample_data.keys():
        plot_content.update({s_name: {}})
        # Read 1
        cycles = sample_data[s_name]["Reads"][0]["Cycles"]
        R1_cycle_num = len(cycles)
        r1r2_split = max(r1r2_split, R1_cycle_num)

    for s_name in sample_data.keys():
        plot_content.update({s_name: {}})
        # Read 1
        cycles = sample_data[s_name]["Reads"][0]["Cycles"]
        R1_cycle_num = len(cycles)
        for cycle in cycles:
            cycle_no = int(cycle["Cycle"])
            adapter_percent = cycle["PercentReadsTrimmed"]
            plot_content[s_name].update({cycle_no: adapter_percent})
        # Read 2
        cycles = sample_data[s_name]["Reads"][1]["Cycles"]
        R2_cycle_num = len(cycles)
        for cycle in cycles:
            cycle_no = int(cycle["Cycle"]) + r1r2_split
            adapter_percent = cycle["PercentReadsTrimmed"]
            plot_content[s_name].update({cycle_no: adapter_percent})
    config = {
        "description": "adapter content",
        "xlab": "Cycle",
        "ylab": "Percentage",
        "xPlotLines": [{"color": "#FF0000", "width": 2, "value": r1r2_split, "dashStyle": "Dash"}],
        "ymax": 100,
        "id": "per_cycle_adapter_content",
        "title": "bases2fastq: per cycle adapter content",
        "ylab": "Percentage",
    }
    plot_name = "Per Sample Adapter Content"
    config.update({"colors": sample_color})
    plot_html = linegraph.plot(plot_content, pconfig=config)
    anchor = "adapter_content"
    description = "Plot of adapter content by cycle"
    helptext = """
    This section plots the adapter content percentage by cycles
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content
