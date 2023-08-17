from multiqc.plots import bargraph, table, linegraph, scatter
import numpy as np
import re
from collections import OrderedDict
from multiqc.utils import report
import json

"""
Functions for plotting per sample information of bases2fastq
"""


def tabulate_sample_stats(sampleData, groupLookupDict, sampleColor):
    """
    Tabulate general information and statistics per sample
    """
    plotContent = dict()
    for s_name in sampleData.keys():
        general_stats = dict()
        general_stats.update({"Run Name": sampleData[s_name]["RunName"]})
        general_stats.update({"Assigned Polonies": sampleData[s_name]["NumPolonies"]})
        general_stats.update({"Sample Percentage Q30": sampleData[s_name]["PercentQ30"]})
        general_stats.update({"Sample Percentage Q40": sampleData[s_name]["PercentQ40"]})
        general_stats.update({"Sample Average Base Quality": sampleData[s_name]["QualityScoreMean"]})
        general_stats.update({"Assigned Yield(Gb)": sampleData[s_name]["Yield"]})
        general_stats.update({"Sample Read1 Average Length": sampleData[s_name]["Reads"][0]["MeanReadLength"]})
        general_stats.update({"Sample Read2 Average Length": sampleData[s_name]["Reads"][1]["MeanReadLength"]})
        general_stats.update({"Group": groupLookupDict[s_name]})
        plotContent.update({s_name: general_stats})
    headers = {header: {"title": header} for header in general_stats.keys()}
    headers["Assigned Polonies"].update({"format": "{d}"})
    config = {"description": "Table of per sample key informations", "no_beeswarm": True}
    plotName = "Sample QC metrics table"
    plotHtml = table.plot(plotContent, headers, pconfig=config)
    anchor = "sample_qc_metrics_table"
    description = "table of general QC metrics by sample"
    helptext = """
    This section shows numbers of some metrics that indicate the quality of each sa,[;e]. 
    
    Polony numbers: total number of polonies assigned to this sample. Each polony yields one read1 and one read2
    
    Percentage Q30: percentage of bases that has >=30 base quality (error rate <= 10^-3) 
    
    Percentage Q40: percentage of bases that has >=40 base quality (error rate <= 10^-4)
    
    Average Base Quality: average base qualities across all bases
    
    Yield: total volume of data that has been assigned to any sample
    
    Read1 Average Length: the average length of Read 1
    
    Read2 Average Length: the average length of Read 2
    
    Group: A group tag for assigning colors in the plot. The default group of each sample will be their sequencing
    run name. To customize group tags, you can either
    
    1) Set the project name when running bases2fastq. In this case the group tags will be project name.
    
    2) Generate a csv file that has the columns "Run Name","Sample Name" and "Group". An easy way of 
    generating this csv will be run multiqc with default grouping, copy this table into a csv file, 
    and then change group name manually. Then move the csv file into any location under the analyzed 
    directory.
    """
    return plotHtml, plotName, anchor, description, helptext


def sequence_content_plot(sampleData, groupLookupDict, colorDict):
    """Create the epic HTML for the FastQC sequence content heatmap"""

    # Prep the data
    data = dict()

    r1r2_split = 0
    for s_name in sorted(sampleData.keys()):
        data[s_name] = {}
        R1 = sampleData[s_name]["Reads"][0]["Cycles"]
        r1r2_split = max(r1r2_split, len(R1))


    for s_name in sorted(sampleData.keys()):
        data[s_name] = {}
        R1 = sampleData[s_name]["Reads"][0]["Cycles"]
        for cycle in range(len(R1)):
            baseNo = str(cycle + 1)
            data[s_name].update({baseNo: {base.lower(): R1[cycle]["BaseComposition"][base] for base in "ACTG"}})
            data[s_name][baseNo]["base"] = baseNo
            tot = sum([data[s_name][baseNo][base] for base in ["a", "c", "t", "g"]])
            for base in ["a", "c", "t", "g"]:
                data[s_name][baseNo][base] = (float(data[s_name][baseNo][base]) / float(tot)) * 100.0

        R2 = sampleData[s_name]["Reads"][1]["Cycles"]
        for cycle in range(len(R2)):
            baseNo = str(cycle + 1 + r1r2_split)
            data[s_name].update({baseNo: {base.lower(): R2[cycle]["BaseComposition"][base] for base in "ACTG"}})
            data[s_name][baseNo]["base"] = baseNo
            tot = sum([data[s_name][baseNo][base] for base in ["a", "c", "t", "g"]])
            for base in ["a", "c", "t", "g"]:
                data[s_name][baseNo][base] = (float(data[s_name][baseNo][base]) / float(tot)) * 100.0
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
        c=json.dumps(["bases2fastq", colorDict]),
        r=r1r2_split,
    )
    plotName = "Per Cycle Sequence Composition"
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
    return html, plotName, anchor, description, helptext


def plot_per_cycle_N_content(sampleData, groupLookupDict, colorDict):
    data = dict()
    r1r2_split = 0 
    for s_name in sorted(sampleData.keys()):
        data[s_name] = {}
        R1 = sampleData[s_name]["Reads"][0]["Cycles"]
        R1CycleNum = len(R1)
        r1r2_split = max(r1r2_split, R1CycleNum)

    for s_name in sorted(sampleData.keys()):
        data[s_name] = {}
        R1 = sampleData[s_name]["Reads"][0]["Cycles"]
        R1CycleNum = len(R1)
        for cycle in range(len(R1)):
            baseNo = str(cycle + 1)
            if sum(R1[cycle]["BaseComposition"].values()) == 0:
                data[s_name].update({baseNo: 0})
            else:
                data[s_name].update(
                    {baseNo: R1[cycle]["BaseComposition"]["N"] / sum(R1[cycle]["BaseComposition"].values()) * 100.0}
                )

        R2 = sampleData[s_name]["Reads"][1]["Cycles"]
        R2CycleNum = len(R2)
        for cycle in range(len(R2)):
            baseNo = str(cycle + 1 + r1r2_split)
            if sum(R2[cycle]["BaseComposition"].values()) == 0:
                data[s_name].update({baseNo: 0})
            else:
                data[s_name].update(
                    {baseNo: R2[cycle]["BaseComposition"]["N"] / sum(R2[cycle]["BaseComposition"].values()) * 100.0}
                )

    plotContent = data
    config = {
        "xlab": "cycle",
        "ylab": "Base Quality",
        "ymax": 100,
        "xPlotLines": [{"color": "#FF0000", "width": 2, "value": r1r2_split, "dashStyle": "Dash"}],
        "colors": colorDict,
        "ymin": 0,
        "id":"per_cycle_n_content",
        "title":"bases2fastq: per cycle N content percentage",
        "ylab":"Percentage"
    }
    plotHtml = linegraph.plot(plotContent, pconfig=config)
    plotName = "Per Cycle N Content"
    anchor = "n_content"
    description = "N content by cycle"
    helptext = """
    This section plots the percentage of unidentified bases ("N" bases) by each sequencing cycle. 
    Read 1 and Read 2 are separated by a red dashed line
    """
    return plotHtml, plotName, anchor, description, helptext


def plot_per_read_gc_hist(sampleData, groupLookupDict, sampleColor):
    """
    Plot gc histogram per sample
    """
    gcHistDict = dict()
    for s_name in sampleData.keys():
        R1GcCounts = sampleData[s_name]["Reads"][0]["PerReadGCCountHistogram"]
        R2GcCounts = sampleData[s_name]["Reads"][1]["PerReadGCCountHistogram"]
        R1R2GcCounts = [r1 + r2 for r1, r2 in zip(R1GcCounts, R2GcCounts)]
        totalReads = sum(R1R2GcCounts)
        gcHistDict.update({s_name: {}})
        RLen = len(R1GcCounts)
        for gc in range(0, RLen):
            gcHistDict[s_name].update({gc / RLen * 100: R1R2GcCounts[gc] / totalReads * 100})

    # perReadQualityHistogram
    plotContent = gcHistDict
    plot_function = linegraph.plot

    config = {"description": "GC", 
            "xlab": "GC content", 
            "ylab": "Percentage", 
            "colors": sampleColor,
            "id":"gc_hist",
            "title":"bases2fastq: per sample GC content histogram",
            "ylab":"Percentage"
    }
    plotName = "per sample GC histogram"
    plotHtml = linegraph.plot(plotContent, pconfig=config)
    anchor = "gc_histogram"
    description = "The histogram of distributions of percentage GC in each read"
    helptext = """
    This section plots the distribution of percentage GC in each reads (range: 0-100)
    """
    return plotHtml, plotName, anchor, description, helptext


def plot_adapter_content(sampleData, groupLookupDict, sampleColor):
    """
    Plot adapter content per sample
    """
    plotContent = dict()
    
    r1r2_split = 0
    for s_name in sampleData.keys():
        plotContent.update({s_name: {}})
        # Read 1
        cycles = sampleData[s_name]["Reads"][0]["Cycles"]
        R1CycleNum = len(cycles)
        r1r2_split = max(r1r2_split, R1CycleNum)

    for s_name in sampleData.keys():
        plotContent.update({s_name: {}})
        # Read 1
        cycles = sampleData[s_name]["Reads"][0]["Cycles"]
        R1CycleNum = len(cycles)
        for cycle in cycles:
            cycleNo = int(cycle["Cycle"])
            adapterPercent = cycle["PercentReadsTrimmed"]
            plotContent[s_name].update({cycleNo: adapterPercent})
        # Read 2
        cycles = sampleData[s_name]["Reads"][1]["Cycles"]
        R2CycleNum = len(cycles)
        for cycle in cycles:
            cycleNo = int(cycle["Cycle"]) + r1r2_split
            adapterPercent = cycle["PercentReadsTrimmed"]
            plotContent[s_name].update({cycleNo: adapterPercent})
    config = {
        "description": "adapter content",
        "xlab": "Cycle",
        "ylab": "Percentage",
        "xPlotLines": [{"color": "#FF0000", "width": 2, "value": r1r2_split, "dashStyle": "Dash"}],
        "ymax": 100,
        "id":"per_cycle_adapter_content",
        "title":"bases2fastq: per cycle adapter content",
        "ylab":"Percentage"
    }
    plotName = "Per Sample Adapter Content"
    config.update({"colors": sampleColor})
    plotHtml = linegraph.plot(plotContent, pconfig=config)
    anchor = "adapter_content"
    description = "The plot of adapter content by cycle"
    helptext = """
    This section plots the adapter content percentage by cycles
    """
    return plotHtml, plotName, anchor, description, helptext
