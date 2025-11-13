from multiqc.plots import linegraph, table
from multiqc import config

"""
Functions for plotting per sample information of bases2fastq
"""


def tabulate_sample_stats(sample_data, group_lookup_dict, project_lookup_dict, sample_color):
    """
    Tabulate general information and statistics per sample
    """
    plot_content = dict()
    for s_name in sample_data.keys():
        general_stats = dict()
        general_stats.update({"group": group_lookup_dict[s_name]})
        general_stats.update({"project": project_lookup_dict.get(s_name, "")})
        general_stats.update({"num_polonies_sample": sample_data[s_name]["NumPolonies"]})
        general_stats.update({"yield_sample": sample_data[s_name]["Yield"]})
        general_stats.update({"mean_base_quality_sample": sample_data[s_name]["QualityScoreMean"]})
        general_stats.update({"percent_q30_sample": sample_data[s_name]["PercentQ30"]})
        general_stats.update({"percent_q40_sample": sample_data[s_name]["PercentQ40"]})
        plot_content.update({s_name: general_stats})

    headers = {}

    headers["group"] = {
        "title": "Group",
        "description": "Run/Sample group label.",
        "bgcols": sample_color,
        "scale": False,
    }
    headers["project"] = {
        "title": "Project",
        "description": "(optional) Run/Sample project label.",
        "bgcols": sample_color,
        "scale": False,
    }
    headers["num_polonies_sample"] = {
        "title": f"# Polonies ({config.base_count_prefix})",
        "description": f"The total number of polonies that are calculated for the run. ({config.base_count_desc})",
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
        "description": "The percentage of ≥ Q30 Q scores for the sample. This includes assigned reads and excludes filtered reads and no calls",
        "max": 100,
        "min": 0,
        "scale": "RdYlGn",
        "suffix": "%",
    }
    headers["percent_q40_sample"] = {
        "title": "Percent Q40",
        "description": "The percentage of ≥ Q40 Q scores for the sample. This includes assigned reads and excludes filtered reads and no calls",
        "max": 100,
        "min": 0,
        "scale": "RdYlGn",
        "suffix": "%",
    }

    pconfig = {"id": "sample_qc_metric_table", "title": "Sample QC Metrics Table", "no_violin": True}

    plot_name = "Sample QC Metrics Table"
    plot_html = table.plot(plot_content, headers, pconfig=pconfig)
    anchor = "sample_qc_metrics_table"
    description = "QC metrics per unique sample"
    helptext = """
     This section displays metrics that indicate the quality of each sample: \n
       - Sample Name: Unique identifier composed of (RunName)__(UUID)__(SampleName), where (RunName) maps to the AVITI run name, (UUID) maps to the unique Bases2Fastq analysis result, and (SampleName) maps to the sample name as specified in the RunManifest.csv.
       - Group: Run/Sample group label that assigns colors in the plot.  To customize group tags:\n
           - 1) Set the project name when running Bases2Fastq. In this case the group tags will be project name.\n
           - 2) Generate a csv file with the suffix "_b2fgroup.csv", containing the columns "Sample Name" and "Group".\n
       - Number of Polonies: The total number of polonies that are assigned to the sample.\n
       - Assigned Yield (Gb): The sample yield that is based on assigned reads in gigabases.\n
       - Quality Score Mean: The average  Q score of base calls for the sample.\n
       - Percent Q30: The percentage of ≥ Q30 Q scores for the sample. This includes assigned reads and excludes filtered reads and no calls.\n
       - Percent Q40: The percentage of ≥ Q40 Q scores for the sample. This includes assigned reads and excludes filtered reads and no calls\n
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content


def sequence_content_plot(sample_data, group_lookup_dict, project_lookup_dict, color_dict):
    """Create the epic HTML for the FastQC sequence content heatmap"""

    # Prep the data
    data = dict()

    r1r2_split = 0
    for s_name in sorted(sample_data.keys()):
        paired_end = True if len(sample_data[s_name]["Reads"]) > 1 else False
        for base in "ACTG":
            base_s_name = "__".join([s_name, base])
            data[base_s_name] = {}
            R1 = sample_data[s_name]["Reads"][0]["Cycles"]
            r1r2_split = max(r1r2_split, len(R1))

    for s_name in sorted(sample_data.keys()):
        R1 = sample_data[s_name]["Reads"][0]["Cycles"]
        for cycle in range(len(R1)):
            base_no = cycle + 1

            tot = sum([R1[cycle]["BaseComposition"][base] for base in ["A", "C", "T", "G"]])

            for base in "ACTG":
                base_s_name = "__".join([s_name, base])
                data[base_s_name].update(
                    {base_no: float(R1[cycle]["BaseComposition"][base] / float(tot)) * 100.0 if tot > 0 else None}
                )

        if paired_end:
            R2 = sample_data[s_name]["Reads"][1]["Cycles"]
            for cycle in range(len(R2)):
                base_no = cycle + 1 + r1r2_split
                tot = sum([R2[cycle]["BaseComposition"][base] for base in ["A", "C", "T", "G"]])

                for base in "ACTG":
                    base_s_name = "__".join([s_name, base])
                    data[base_s_name].update(
                        {base_no: float(R2[cycle]["BaseComposition"][base] / float(tot)) * 100.0 if tot > 0 else None}
                    )

    plot_content = data

    pconfig = {
        "xlab": "cycle",
        "ylab": "Percentage",
        "x_lines": [{"color": "#FF0000", "width": 2, "value": r1r2_split, "dashStyle": "dash"}],
        "colors": color_dict,
        "ymin": 0,
        "id": "per_cycle_base_content",
        "title": "bases2fastq: Per Cycle Base Content Percentage",
    }
    plot_html = linegraph.plot(plot_content, pconfig=pconfig)
    plot_name = "Per Cycle Base Content"
    anchor = "base_content"
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


def plot_per_cycle_N_content(sample_data, group_lookup_dict, project_lookup_dict, color_dict):
    data = dict()
    r1r2_split = 0
    for s_name in sorted(sample_data.keys()):
        data[s_name] = {}
        R1 = sample_data[s_name]["Reads"][0]["Cycles"]
        R1_cycle_num = len(R1)
        r1r2_split = max(r1r2_split, R1_cycle_num)

    for s_name in sorted(sample_data.keys()):
        paired_end = True if len(sample_data[s_name]["Reads"]) > 1 else False
        R1 = sample_data[s_name]["Reads"][0]["Cycles"]
        R1_cycle_num = len(R1)
        for cycle in range(len(R1)):
            base_no = cycle + 1
            if sum(R1[cycle]["BaseComposition"].values()) == 0:
                data[s_name].update({base_no: 0})
            else:
                data[s_name].update(
                    {base_no: R1[cycle]["BaseComposition"]["N"] / sum(R1[cycle]["BaseComposition"].values()) * 100.0}
                )

        if paired_end:
            R2 = sample_data[s_name]["Reads"][1]["Cycles"]
            for cycle in range(len(R2)):
                base_no = cycle + 1 + r1r2_split
                if sum(R2[cycle]["BaseComposition"].values()) == 0:
                    data[s_name].update({base_no: 0})
                else:
                    data[s_name].update(
                        {
                            base_no: R2[cycle]["BaseComposition"]["N"]
                            / sum(R2[cycle]["BaseComposition"].values())
                            * 100.0
                        }
                    )

    plot_content = data
    pconfig = {
        "xlab": "cycle",
        "ylab": "Percentage",
        "x_lines": [{"color": "#FF0000", "width": 2, "value": r1r2_split, "dashStyle": "dash"}],
        "colors": color_dict,
        "ymin": 0,
        "ymax": 100,
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


def plot_per_read_gc_hist(sample_data, group_lookup_dict, project_lookup_dict, sample_color):
    """
    Plot GC Histogram per Sample
    """
    gc_hist_dict = dict()
    for s_name in sample_data.keys():
        R1_gc_counts = sample_data[s_name]["Reads"][0]["PerReadGCCountHistogram"]
        R2_gc_counts = [0] * len(R1_gc_counts)
        if len(sample_data[s_name]["Reads"]) > 1:
            R2_gc_counts_raw = sample_data[s_name]["Reads"][1]["PerReadGCCountHistogram"]
            # Handle potential length mismatch between R1 and R2 GC counts
            if len(R2_gc_counts_raw) == len(R1_gc_counts):
                R2_gc_counts = R2_gc_counts_raw
            else:
                # Pad shorter list with zeros or truncate longer list
                min_len = min(len(R1_gc_counts), len(R2_gc_counts_raw))
                R2_gc_counts = R2_gc_counts_raw[:min_len] + [0] * (len(R1_gc_counts) - min_len)

        R1R2_gc_counts = [r1 + r2 for r1, r2 in zip(R1_gc_counts, R2_gc_counts)]
        totalReads = sum(R1R2_gc_counts)
        gc_hist_dict.update({s_name: {}})

        if totalReads > 0 and len(R1R2_gc_counts) > 0:  # Avoid division by zero and empty data
            RLen = len(R1_gc_counts)
            for gc in range(len(R1R2_gc_counts)):
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


def plot_adapter_content(sample_data, group_lookup_dict, project_lookup_dict, sample_color):
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
        paired_end = True if len(sample_data[s_name]["Reads"]) > 1 else False
        plot_content.update({s_name: {}})
        # Read 1
        cycles = sample_data[s_name]["Reads"][0]["Cycles"]
        for cycle in cycles:
            cycle_no = int(cycle["Cycle"])
            adapter_percent = cycle["PercentReadsTrimmed"]
            plot_content[s_name].update({cycle_no: adapter_percent})
        # Read 2
        if paired_end:
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
