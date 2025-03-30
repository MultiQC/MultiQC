import json
from collections import OrderedDict

from multiqc.plots import linegraph, table
from multiqc.utils import config, report

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
        general_stats.update({"group": group_lookup_dict[s_name]})
        general_stats.update({"num_polonies_sample": sample_data[s_name]["NumPolonies"]})
        general_stats.update({"yield_sample": sample_data[s_name]["Yield"]})
        general_stats.update({"mean_base_quality_sample": sample_data[s_name]["QualityScoreMean"]})
        general_stats.update({"percent_q30_sample": sample_data[s_name]["PercentQ30"]})
        general_stats.update({"percent_q40_sample": sample_data[s_name]["PercentQ40"]})
        plot_content.update({s_name: general_stats})

    headers = OrderedDict()

    headers["group"] = {
        "title": "Group",
        "description": "Run/Sample group label.",
        "bgcols": sample_color,
        "scale": False,
    }
    headers["num_polonies_sample"] = {
        "title": f"# Polonies ({config.base_count_prefix})",
        "description": f"The (total) number of polonies calculated for the sample ({config.base_count_desc})",
        "min": 0,
        "scale": "Blues",
        "shared_key": "base_count",
    }
    headers["yield_sample"] = {
        "title": "Yield (Gb)",
        "description": "The sample yield based on assigned reads in gigabases",
        "scale": "Greens",
    }
    headers["mean_base_quality_sample"] = {
        "title": "Mean Base Quality",
        "description": "Average base quality across R1/R2",
        "min": 0,
        "scale": "Spectral",
    }
    headers["percent_q30_sample"] = {
        "title": "Percent Q30",
        "description": "The percentage of ≥ Q30 (base) Q-scores for the run, including assigned and unassigned reads",
        "max": 100,
        "min": 0,
        "scale": "RdYlGn",
        "suffix": "%",
    }
    headers["percent_q40_sample"] = {
        "title": "Percent Q40",
        "description": "The percentage of ≥ Q40 (base) Q-scores for the run, including assigned and unassigned reads",
        "max": 100,
        "min": 0,
        "scale": "RdYlGn",
        "suffix": "%",
    }

    pconfig = {
            "id": "sample_qc_metric_table",
            "title": "Sample QC Metrics Table",
            "no_violin": True
        }

    plot_name = "Sample QC metrics table"
    plot_html = table.plot(plot_content, headers, pconfig=pconfig)
    anchor = "sample_qc_metrics_table"
    description = "table of general QC metrics by sample"
    helptext = """
    This section shows numbers of some metrics that indicate the quality of each sequencing run: \n
       - Sample Name: Name showing the (RunName)__(UUID)__(SampleName).  (RunName) maps to the AVITI run name.  (UUID) maps to the unique bases2fastq analysis result.  (SampleName) maps to the sample name as specified in the RunManifest.csv.
       - Group: Run/Sample group label for assigning colors in the plot.  To customize group tags, you can:\n
           - 1) Set the project name when running bases2fastq. In this case the group tags will be project name.\n
           - 2) Generate a csv file with the suffix "_b2fgroup.csv", containing the columns "Sample Name" and "Group".\n
       - Number of Polonies: The (total) number of polonies calculated for the run\n
       - Assigned Yield (Gb): The run yield based on assigned reads in gigabases\n
       - Quality Score Mean: The average Q-score of base calls for a sample\n
       - Percent Q30: The percentage of ≥ Q30 (base) Q-scores for the run, including assigned and unassigned reads\n
       - Percent Q40: The percentage of ≥ Q40 (base) Q-scores for the run, including assigned and unassigned reads
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
    The heatmap is formatted in a similar way as FastQC per base sequence content. The differences are:
    1) The pass/fail flag in the original FastQC module are replaced with a color reflecting the group tag of the samples
    2) Read1 and Read2 is co-plotted and separate by a white band in the heatmap and a red dashed in each individual sample plot.

    The following is a copy of helptext of the FastQC module:
    
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
        for cycle in range(len(R2)):
            base_no = str(cycle + 1 + r1r2_split)
            if sum(R2[cycle]["BaseComposition"].values()) == 0:
                data[s_name].update({base_no: 0})
            else:
                data[s_name].update(
                    {base_no: R2[cycle]["BaseComposition"]["N"] / sum(R2[cycle]["BaseComposition"].values()) * 100.0}
                )

    plot_content = data
    pconfig = {
        "xlab": "cycle",
        "ylab": "Percentage",
        "ymax": 100,
        "x_lines": [{"color": "#FF0000", "width": 2, "value": r1r2_split, "dashStyle": "dash"}],
        "colors": color_dict,
        "ymin": 0,
        "id": "per_cycle_n_content",
        "title": "bases2fastq: Per Cycle N Content Percentage",
    }
    plot_html = linegraph.plot(plot_content, pconfig=pconfig)
    plot_name = "Per Cycle N Content"
    anchor = "n_content"
    description = """
    Percentage of unidentified bases ("N" bases) by each sequencing cycle.
    Read 1 and Read 2 are separated by a red dashed line
    """
    helptext = """
    If a sequencer is unable to make a base call with sufficient confidence then it will
    normally substitute an `N` rather than a conventional base call. This graph shows the
    percentage of base calls at each position for which an `N` was called.

    It's not unusual to see a very low proportion of Ns appearing in a sequence, especially
    nearer the end of a sequence. However, if this proportion rises above a few percent
    it suggests that the analysis pipeline was unable to interpret the data well enough to
    make valid base calls.
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content


def plot_per_read_gc_hist(sample_data, group_lookup_dict, sample_color):
    """
    Plot GC Histogram per Sample
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

    pconfig = {
        "xlab": "% GC",
        "ylab": "Percentage",
        "colors": sample_color,
        "id": "gc_hist",
        "title": "bases2fastq: Per Sample GC Content Histogram",
    }
    plot_name = "Per Sample GC Histogram"
    plot_html = linegraph.plot(plot_content, pconfig=pconfig)
    anchor = "gc_histogram"
    description = "Histogram of distributions of percentage GC in each read"
    helptext = """
    GC content across the whole length of each sequence.
    
    In a normal random library you would expect to see a roughly normal distribution
    of GC content where the central peak corresponds to the overall GC content of
    the underlying genome. Since we don't know the GC content of the genome the
    modal GC content is calculated from the observed data and used to build a
    reference distribution.

    An unusually shaped distribution could indicate a contaminated library or
    some other kinds of biased subset. A normal distribution which is shifted
    indicates some systematic bias which is independent of base position. If there
    is a systematic bias which creates a shifted normal distribution then this won't
    be flagged as an error by the module since it doesn't know what your genome's
    GC content should be.
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content


def plot_adapter_content(sample_data, group_lookup_dict, sample_color):
    """
    Plot Adapter Content per Sample
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
        for cycle in cycles:
            cycle_no = int(cycle["Cycle"])
            adapter_percent = cycle["PercentReadsTrimmed"]
            plot_content[s_name].update({cycle_no: adapter_percent})
        # Read 2
        cycles = sample_data[s_name]["Reads"][1]["Cycles"]
        for cycle in cycles:
            cycle_no = int(cycle["Cycle"]) + r1r2_split
            adapter_percent = cycle["PercentReadsTrimmed"]
            plot_content[s_name].update({cycle_no: adapter_percent})
    pconfig = {
        "id": "per_cycle_adapter_content",
        "title": "bases2fastq: Per Cycle Adapter Content",
        "xlab": "Cycle",
        "ylab": "% of Sequences",
        "x_lines": [{"color": "#FF0000", "width": 2, "value": r1r2_split, "dashStyle": "dash"}],
        "ymax": 100,
    }
    plot_name = "Per Sample Adapter Content"
    pconfig.update({"colors": sample_color})
    plot_html = linegraph.plot(plot_content, pconfig=pconfig)
    anchor = "adapter_content"
    description = "Adapter content per cycle"
    helptext = """
    The plot shows a cumulative percentage count of the proportion
    of your library which has seen each of the adapter sequences at each cycle.
    Once a sequence has been seen in a read, it is counted as being present
    right through to the end of the read, so the percentages you see will only
    increase as the read length goes on.
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content
