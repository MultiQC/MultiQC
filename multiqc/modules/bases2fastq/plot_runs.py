import math

from multiqc.plots import bargraph, linegraph, table
from multiqc import config


"""
Functions for plotting per run information of bases2fastq
"""


def plot_run_stats(run_data, color_dict):
    """
    Plot a bar graph for polony numbers, Q30/Q40, index assignment rate and yields for each run
    """
    run_names = list(run_data.keys())
    run_names.sort()
    num_polonies = dict()
    yields = dict()
    for run in run_names:
        # Index Assignment Polonies and Yields ###
        # percent_assigned = run_data[run].get("PercentAssignedReads",100.0)
        percent_assigned = run_data[run]["PercentAssignedReads"]

        percent_perfect_assigned = (
            100.00 - run_data[run]["PercentMismatch"]
        )  # percentage of assigned polonies that has perfect index pairs
        percent_perfect_total = (
            0.01 * percent_assigned * percent_perfect_assigned
        )  # percentage of total polonies that has perfect index pairs
        percent_imperfect_total = (
            percent_assigned - percent_perfect_total
        )  # percentage of total polonies that is assigned but has index mismatch

        num_polonies_run = {}
        polonies = run_data[run]["NumPolonies"]
        num_polonies_run["Perfect Index"] = math.ceil(polonies * percent_perfect_total * 0.01)
        num_polonies_run["Mismatched Index"] = math.ceil(polonies * percent_imperfect_total * 0.01)
        num_polonies_run["Unassigned"] = (
            polonies - num_polonies_run["Perfect Index"] - num_polonies_run["Mismatched Index"]
        )
        num_polonies[run] = num_polonies_run

        total_yield_run = {}
        total_yield = run_data[run].get("TotalYield", 300.0)
        total_yield_run["Perfect Index"] = total_yield * percent_perfect_total * 0.01
        total_yield_run["Mismatched Index"] = total_yield * percent_imperfect_total * 0.01
        total_yield_run["Unassigned"] = (
            total_yield - total_yield_run["Perfect Index"] - total_yield_run["Mismatched Index"]
        )
        yields[run] = total_yield_run

    plot_content = [num_polonies, yields]
    pconfig = {
        "data_labels": [
            {"name": "Number of Polonies", "ylab": "Number of Polonies", "format": "{d}"},
            {"name": "Yield (Gb)", "ylab": "Gb"},
        ],
        "cpswitch": True,
        "stacking": "normal",
        "id": "run_metrics_bar",
        "title": "bases2fastq: General Sequencing Run QC metrics plot",
        "ylab": "QC",
    }
    cats = [
        {
            "Perfect Index": {"name": "Perfect Index", "color": "#7cb4ec"},
            "Mismatched Index": {"name": "Mismatched Index", "color": "#ff7518"},
            "Unassigned": {"name": "Unassigned Index", "color": "#434348"},
        }
    ] * 2

    plot_name = "Sequencing Run Yield"
    plot_html = bargraph.plot(plot_content, cats, pconfig=pconfig)
    anchor = "run_yield_plot"
    description = "Bar plots of sequencing run yields. Please see individual run reports for details"
    helptext = """
    This section shows and compare the yield and index assignment rate of each sequencing run.\n\n
       - Number of Polonies: The total number of polonies that are calculated for the run.\n
       - Yield: The total yield of all assigned reads in gigabases.
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content


def tabulate_run_stats(run_data, color_dict):
    """
    Tabulate general information and statistics of each run
    """
    plot_content = dict()
    for s_name in run_data.keys():
        run_stats = dict()
        run_stats.update({"num_polonies_run": int(run_data[s_name]["NumPolonies"])})
        run_stats.update({"percent_assigned_run": run_data[s_name].get("PercentAssignedReads", 100.0)})
        run_stats.update({"yield_run": run_data[s_name]["AssignedYield"]})
        run_stats.update({"mean_base_quality_run": run_data[s_name]["QualityScoreMean"]})
        run_stats.update({"percent_q30_run": run_data[s_name]["PercentQ30"]})
        run_stats.update({"percent_q40_run": run_data[s_name]["PercentQ40"]})
        plot_content.update({s_name: run_stats})

    headers = {}
    headers["num_polonies_run"] = {
        "title": f"# Polonies ({config.base_count_prefix})",
        "description": f"The total number of polonies that are calculated for the run. ({config.base_count_desc})",
        "min": 0,
        "scale": "RdYlGn",
        "shared_key": "base_count",
    }
    headers["percent_assigned_run"] = {
        "title": "% Assigned Reads",
        "description": "The percentage of reads assigned to sample(s)",
        "max": 100,
        "min": 0,
        "scale": "BuPu",
        "suffix": "%",
    }
    headers["yield_run"] = {
        "title": "Yield (Gb)",
        "description": "The run yield based on assigned reads in gigabases",
        "scale": "Blues",
    }
    headers["mean_base_quality_run"] = {
        "title": "Quality Score Mean",
        "description": "The average base quality across Read 1 and Read 2.",
        "min": 0,
        "scale": "Spectral",
    }
    headers["percent_q30_run"] = {
        "title": "Percent Q30",
        "description": "The percentage of ≥ Q30 Q scores for the run. This includes assigned and unassigned reads and excludes filtered reads and no calls.",
        "max": 100,
        "min": 0,
        "scale": "RdYlGn",
        "suffix": "%",
    }
    headers["percent_q40_run"] = {
        "title": "Percent Q40",
        "description": "The percentage of ≥ Q40 Q scores for the run. This includes assigned and unassigned reads and excludes filtered reads and no calls.",
        "max": 100,
        "min": 0,
        "scale": "RdYlGn",
        "suffix": "%",
    }

    pconfig = {
        "title": "bases2fastq: General Sequencing Run QC metrics",
        "col1_header": "Run Name",
        "id": "run_metrics_table",
        "ylab": "QC",
    }

    plot_name = "Sequencing Run QC Metrics Table"
    plot_html = table.plot(plot_content, headers, pconfig=pconfig)
    anchor = "run_qc_metrics_table"
    description = "QC metrics per run"
    helptext = """
    This section displays metrics that indicate the quality of each sequencing run: \n
       - Run Name: Unique identifier composed of (RunName)__(UUID), where (RunName) maps to the AVITI run name and (UUID) maps to the unique Bases2Fastq analysis result.\n
       - Number of Polonies: The total number of polonies that are calculated for the run.\n
       - Percentage Assigned Reads: The percentage of reads that are assigned to a sample.\n
       - Assigned Yield (Gb): The run yield that is based on assigned reads in gigabases.\n
       - Quality Score Mean: The mean Q score of base calls for the samples. This excludes filtered reads and no calls.\n
       - Percent Q30: The percentage of ≥ Q30 Q scores for the run. This includes assigned and unassigned reads and excludes filtered reads and no calls.\n
       - Percent Q40: The percentage of ≥ Q40 Q scores for the run. This includes assigned and unassigned reads and excludes filtered reads and no calls.\n
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content


def plot_base_quality_hist(run_data, color_dict):
    # Prepare plot data for per base BQ histogram
    bq_hist_dict = dict()
    for s_name in run_data.keys():
        paired_end = True if len(run_data[s_name]["Reads"]) > 1 else False
        R1_base_quality_counts = run_data[s_name]["Reads"][0]["QualityScoreHistogram"]
        R2_base_quality_counts = [0] * len(R1_base_quality_counts)
        if paired_end:
            R2_base_quality_counts = run_data[s_name]["Reads"][1]["QualityScoreHistogram"]
        R1R2_base_quality_counts = [r1 + r2 for r1, r2 in zip(R1_base_quality_counts, R2_base_quality_counts)]
        total_bases = sum(R1R2_base_quality_counts)
        bq_hist_dict.update({s_name: {}})
        for quality in range(0, len(R1R2_base_quality_counts)):
            bq_hist_dict[s_name].update({quality: R1R2_base_quality_counts[quality] / total_bases * 100})

    # Prepare plot data for per read average BQ histogram
    per_read_quality_hist_dict = dict()
    for s_name in run_data.keys():
        paired_end = True if len(run_data[s_name]["Reads"]) > 1 else False
        R1_quality_counts = run_data[s_name]["Reads"][0]["PerReadMeanQualityScoreHistogram"]
        R2_quality_counts = [0] * len(R1_quality_counts)
        if paired_end:
            R2_quality_counts = run_data[s_name]["Reads"][1]["PerReadMeanQualityScoreHistogram"]
        multiplier = 2 if paired_end else 1
        total_reads = run_data[s_name]["NumPolonies"] * multiplier
        R1R2_quality_counts = [r1 + r2 for r1, r2 in zip(R1_quality_counts, R2_quality_counts)]
        per_read_quality_hist_dict.update({s_name: {}})
        for meanQuality in range(0, len(R1R2_quality_counts)):
            per_read_quality_hist_dict[s_name].update(
                {meanQuality: R1R2_quality_counts[meanQuality] / total_reads * 100}
            )
    plot_content = [bq_hist_dict, per_read_quality_hist_dict]

    # Config for switching dataset
    pconfig = {
        "data_labels": [
            {
                "name": "Quality Per Base",
                "description": "Histogram of bases quality",
                "ymin": 0,
                "ylabel": "Percentage of base quality",
                "xlabel": "base quality",
                "colors": color_dict,
            },
            {
                "name": "Qualiter Per Read",
                "description": "Histogram of average read base quality",
                "ymin": 0,
                "ylabel": "Percentage of read quality",
                "xlabel": "base quality",
                "colors": color_dict,
            },
        ],
        "id": "per_run_bq_hist",
        "title": "bases2fastq: Quality Histograms",
        "ylab": "Percentage",
    }
    plot_html = linegraph.plot(plot_content, pconfig=pconfig)
    plot_name = "Run Base Quality Histogram"
    anchor = "bq_hist"
    description = "Histogram of run base qualities"
    helptext = """
    Run base qualities histogram, summarised by bases and reads. 
    Use tabs to switch between the views:\n
       - Quality Per Base: distribution of base qualities.\n
       - Quality Per Read: distribution of read qualities.\n
    \n
    _The y-axis on the graph shows the quality scores. The higher the score, the better
    the base call. The background of the graph divides the y-axis into very good quality
    calls (green), calls of reasonable quality (orange), and calls of poor quality (red).
    The quality of calls on most platforms will degrade as the run progresses, so it is
    common to see base calls falling into the orange area towards the end of a read._
    """

    return plot_html, plot_name, anchor, description, helptext, plot_content


def plot_base_quality_by_cycle(run_data, color_dict):
    # Prepare plot data for median BQ of each cycle

    r1r2_split = 0
    for s_name in run_data.keys():
        paired_end = True if len(run_data[s_name]["Reads"]) > 1 else False
        cycle_dict = dict()
        R1CycleNum = len(run_data[s_name]["Reads"][0]["Cycles"])
        r1r2_split = max(r1r2_split, R1CycleNum)

    median_dict = {}
    for s_name in run_data.keys():
        paired_end = True if len(run_data[s_name]["Reads"]) > 1 else False
        cycle_dict = dict()
        R1CycleNum = len(run_data[s_name]["Reads"][0]["Cycles"])
        for cycle in run_data[s_name]["Reads"][0]["Cycles"]:
            cycle_no = int(cycle["Cycle"])
            cycle_dict.update({cycle_no: cycle["QualityScore50thPercentile"]})
        if paired_end:
            for cycle in run_data[s_name]["Reads"][1]["Cycles"]:
                cycle_no = int(cycle["Cycle"]) + r1r2_split
                cycle_dict.update({cycle_no: cycle["QualityScore50thPercentile"]})
        median_dict.update({s_name: cycle_dict})

    # Prepare plot data for mean BQ of each cycle
    mean_dict = {}
    for s_name in run_data.keys():
        paired_end = True if len(run_data[s_name]["Reads"]) > 1 else False
        # Update each sample cycle info
        cycle_dict = dict()
        for cycle in run_data[s_name]["Reads"][0]["Cycles"]:
            cycle_no = int(cycle["Cycle"])
            cycle_dict.update({cycle_no: cycle["QualityScoreMean"]})
        if paired_end:
            for cycle in run_data[s_name]["Reads"][1]["Cycles"]:
                cycle_no = int(cycle["Cycle"]) + r1r2_split
                cycle_dict.update({cycle_no: cycle["QualityScoreMean"]})
        mean_dict.update({s_name: cycle_dict})

    # Prepare plot data for %Q30 of each cycle
    Q30_dict = {}
    for s_name in run_data.keys():
        paired_end = True if len(run_data[s_name]["Reads"]) > 1 else False
        # Update each sample cycle info
        cycle_dict = dict()
        for cycle in run_data[s_name]["Reads"][0]["Cycles"]:
            cycle_no = int(cycle["Cycle"])
            cycle_dict.update({cycle_no: cycle["PercentQ30"]})
        if paired_end:
            for cycle in run_data[s_name]["Reads"][1]["Cycles"]:
                cycle_no = int(cycle["Cycle"]) + r1r2_split
                cycle_dict.update({cycle_no: cycle["PercentQ30"]})
        Q30_dict.update({s_name: cycle_dict})

    # Prepare plot data for %Q40 of each cycle
    Q40_dict = {}
    for s_name in run_data.keys():
        paired_end = True if len(run_data[s_name]["Reads"]) > 1 else False
        cycle_dict = dict()
        for cycle in run_data[s_name]["Reads"][0]["Cycles"]:
            cycle_no = int(cycle["Cycle"])
            cycle_dict.update({cycle_no: cycle["PercentQ40"]})
        if paired_end:
            for cycle in run_data[s_name]["Reads"][1]["Cycles"]:
                cycle_no = int(cycle["Cycle"]) + r1r2_split
                cycle_dict.update({cycle_no: cycle["PercentQ40"]})
        Q40_dict.update({s_name: cycle_dict})

    # Prepare plot data for % base calls below PF threshold
    below_pf_dict = {}
    for s_name in run_data.keys():
        paired_end = True if len(run_data[s_name]["Reads"]) > 1 else False
        cycle_dict = dict()
        R1CycleNum = len(run_data[s_name]["Reads"][0]["Cycles"])
        if "PercentBelowFilterThreshold" not in run_data[s_name]["Reads"][0]["Cycles"][0]:
            continue
        for cycle in run_data[s_name]["Reads"][0]["Cycles"]:
            cycle_no = int(cycle["Cycle"])
            cycle_dict.update({cycle_no: cycle["PercentBelowFilterThreshold"]})
        if paired_end:
            for cycle in run_data[s_name]["Reads"][1]["Cycles"]:
                cycle_no = int(cycle["Cycle"]) + r1r2_split
                cycle_dict.update({cycle_no: cycle["PercentBelowFilterThreshold"]})
        below_pf_dict.update({s_name: cycle_dict})

    # aggregate plot data
    plot_content = [median_dict, mean_dict, Q30_dict, Q40_dict, below_pf_dict]
    pconfig = {
        "data_labels": [
            {"name": "Median Quality", "xlab": "cycle", "ylab": "Quality"},
            {"name": "Mean Quality", "ylab": "Quality"},
            {"name": "%Q30", "xlab": "cycle", "ylab": "Percentage", "ymax": 100},
            {"name": "%Q40", "xlab": "cycle", "ylab": "Percentage", "ymax": 100},
            {"name": "%Base Calls Below PF", "xlab": "cycle", "ylab": "Percentage", "ymax": 100},
        ],
        "x_lines": [{"color": "#FF0000", "width": 2, "value": r1r2_split, "dashStyle": "dash"}],
        "colors": color_dict,
        "ymin": 0,
        "id": "per_run_quality_by_cycle",
        "title": "bases2fastq: Quality by cycles",
        "ylab": "QC",
    }
    plot_html = linegraph.plot(plot_content, pconfig=pconfig)
    plot_name = "Quality Metrics By Cycle"
    anchor = "per_cycle_quality"
    description = "Per run base qualities by cycle"
    helptext = """
    This section plots the base qualities by each instrument cycle.\n
    Choose between Median Quality, Mean Quality, Percent Q30 or Percentage Q40 per cycle.\n
    Read 1 and Read 2 are separated by a red dashed line.
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content
