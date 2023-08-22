import math
from collections import OrderedDict

import numpy as np

from multiqc.plots import bargraph, linegraph, scatter, table

"""
Functions for plotting per run information of bases2fastq
"""


def plot_run_stats(run_data, color_dict):
    """
    plot bargraph for polony numbers, Q30/Q40, index assignment rate and yields for each run
    """
    run_names = list(run_data.keys())
    run_names.sort()
    num_polonies = dict()
    percent_qualities = dict()
    percent_assigned_dict = dict()
    yields = dict()
    for run in run_names:
        ### Index Assignment Polonies and Yields ###
        percent_assigned = run_data[run]["PercentAssignedReads"]
        percent_unassigned = run_data[run]["PercentMismatch"]
        percent_perfect_assigned = (
            100.00 - run_data[run]["PercentMismatch"]
        )  # percentage of assigned polonies that has perfect index pairs
        percent_perfect_total = (
            0.01 * percent_assigned * percent_perfect_assigned
        )  # percentage of total polonies that has perfect index pairs
        percent_imperfect_total = (
            percent_assigned - percent_perfect_total
        )  # percentage of total polonies that is assigned but has index mismatch

        num_polonies_run = OrderedDict()
        polonies = run_data[run]["NumPolonies"]
        num_polonies_run["Perfect Index"] = math.ceil(polonies * percent_perfect_total * 0.01)
        num_polonies_run["Mismatched Index"] = math.ceil(polonies * percent_imperfect_total * 0.01)
        num_polonies_run["Unassigned"] = (
            polonies - num_polonies_run["Perfect Index"] - num_polonies_run["Mismatched Index"]
        )
        num_polonies[run] = num_polonies_run

        total_yield_run = OrderedDict()
        total_yield = run_data[run]["TotalYield"]
        total_yield_run["Perfect Index"] = total_yield * percent_perfect_total * 0.01
        total_yield_run["Mismatched Index"] = total_yield * percent_imperfect_total * 0.01
        total_yield_run["Unassigned"] = (
            total_yield - total_yield_run["Perfect Index"] - total_yield_run["Mismatched Index"]
        )
        yields[run] = total_yield_run

    plot_content = [num_polonies, yields]
    config = {
        "data_labels": [
            {"name": "Polony Numbers", "ylab": "Number of polonies", "format": "{d}"},
            {"name": "Data Yield", "ylab": "Gb"},
        ],
        "cpswitch": True,
        "stacking": "normal",
        "description": "Bar plots to compare yields of all found sequencing runs.",
        "id": "run_stats_bar",
        "title": "bases2fastq: General Sequencing Run QC stats plot",
        "ylab": "QC",
    }
    cats = [
        {
            "Perfect Index": {"name": "Perfect Index", "color": "#7cb4ec"},
            "Mismatched Index": {"name": "Mismatched Index", "color": "#90ed7d"},
            "Unassigned": {"name": "Unassigned Index", "color": "#434348"},
        }
    ] * 2

    plot_name = "Sequencing Run Yield"
    plot_html = bargraph.plot(plot_content, cats, pconfig=config)
    anchor = "run_yield_plot"
    description = "Bar plots of sequencing run yields. Please see individual run reports for details"
    helptext = """
    This section shows and compare the yield and index assignment rate of each sequencing run. 
    
    Polony numbers: total number of polonies. Each polony will yield one read1 and one read2
    Data Yield: the total volume of data that has been assigned to any sample
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content


def tabulate_run_stats(run_data, color_dict):
    """
    Tabulate general information and statistics of each run
    """
    plot_content = dict()
    for s_name in run_data.keys():
        run_stats = dict()
        run_stats.update({"num_polonies_run": run_data[s_name]["NumPolonies"]})
        run_stats.update({"percent_q30_run": run_data[s_name]["PercentQ30"]})
        run_stats.update({"percent_q40_run": run_data[s_name]["PercentQ40"]})
        run_stats.update({"mean_base_quality": run_data[s_name]["QualityScoreMean"]})
        run_stats.update({"yield_run": run_data[s_name]["AssignedYield"]})
        plot_content.update({s_name: run_stats})

    headers = OrderedDict()

    headers["num_polonies_run"] = {
        "title": "#Polonies",
        "format": "{d}",
        "description": "Percent of reads with perfect index (0 mismatches)",
        "min": 0,
        "scale": "RdYlGn",
    }
    headers["yield_run"] = {"title": "Yield(Gb)", "description": "Total Yield(GB) of run", "scale": "Blues"}
    headers["mean_base_quality_sample"] = {
        "title": "Mean Base Quality",
        "description": "Average base quality across R1/R2",
        "min": 0,
        "scale": "Spectral",
        "suffix": "%",
    }
    headers["percent_q30_run"] = {
        "title": "% Bases Q30",
        "description": "Percent of bases >Q30. Q30 indicates an incorrect base call probability of 1 in 1,000, which equals a 99.9% confidence rate.",
        "max": 100,
        "min": 0,
        "scale": "RdYlGn",
        "suffix": "%",
    }
    headers["percent_q40_run"] = {
        "title": "% Bases Q40",
        "description": "Percent of bases >Q40. Q40 indicates an incorrect base call probability of 1 in 10,000, which equals a 99.99% confidence rate.",
        "max": 100,
        "min": 0,
        "scale": "RdYlGn",
        "suffix": "%",
    }
    headers["percent_assigned"] = {
        "title": "% Bases Assigned",
        "description": "Percent of reads assigned.",
        "max": 100,
        "min": 0,
        "scale": "BuPu",
    }

    config = {
        "title": "bases2fastq: General Sequencing Run QC stats table",
        "descriptions": "Table of per sequencing run key informations",
        "col1_header": "Run Name",
        "id": "run_stats_table",
        "ylab": "QC",
    }

    plot_name = "Sequencing Run QC metrics table"
    plot_html = table.plot(plot_content, headers, pconfig=config)
    anchor = "run_qc_metrics_table"
    description = "Table of general QC metrics"
    helptext = """
    This section shows numbers of some metrics that indicate the quality of each sequencing run. 
    \n
    Polony numbers: total number of polonies. Each polony yields one read1 and one read2
    \n
    Yield: total volume of data that has been assigned to any sample
    \n
    Percentage Q30: percentage of bases that has >=30 base quality (error rate <= 10^-3) 
    \n
    Percentage Q40: percentage of bases that has >=40 base quality (error rate <= 10^-4)
    \n
    Average Base Quality: average base qualities across all bases
    \n
    Percentage assigned: percentage of reads that has been assigned to any sample
    \n
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content


def plot_base_quality_hist(run_data, color_dict):
    # Prepare plot data for per base BQ histogram
    bq_hist_dict = dict()
    for s_name in run_data.keys():
        R1_base_quality_counts = run_data[s_name]["Reads"][0]["QualityScoreHistogram"]
        R2_base_quality_counts = run_data[s_name]["Reads"][1]["QualityScoreHistogram"]
        R1R2_base_quality_counts = [r1 + r2 for r1, r2 in zip(R1_base_quality_counts, R2_base_quality_counts)]
        total_bases = sum(R1R2_base_quality_counts)
        bq_hist_dict.update({s_name: {}})
        for quality in range(0, len(R1R2_base_quality_counts)):
            bq_hist_dict[s_name].update({quality: R1R2_base_quality_counts[quality] / total_bases * 100})

    # Prepare plot data for per read average BQ histogram
    per_read_quality_hist_dict = dict()
    for s_name in run_data.keys():
        R1_quality_counts = run_data[s_name]["Reads"][0]["PerReadMeanQualityScoreHistogram"]
        R2_quality_counts = run_data[s_name]["Reads"][1]["PerReadMeanQualityScoreHistogram"]
        total_reads = run_data[s_name]["NumPolonies"] * 2
        R1R2_quality_counts = [r1 + r2 for r1, r2 in zip(R1_quality_counts, R2_quality_counts)]
        per_read_quality_hist_dict.update({s_name: {}})
        for meanQuality in range(0, len(R1R2_quality_counts)):
            per_read_quality_hist_dict[s_name].update(
                {meanQuality: R1R2_quality_counts[meanQuality] / total_reads * 100}
            )
    plot_content = [bq_hist_dict, per_read_quality_hist_dict]

    # Config for switching dataset
    config = {
        "data_labels": [
            {
                "name": "per base quality histogram",
                "descriptions": "Histogram of bases qualities",
                "ymin": 0,
                "ylabel": "Percentage at base quality",
                "xlabel": "base quality",
                "colors": color_dict,
            },
            {
                "name": "per read quality histogram",
                "descriptions": "Histogram of read average base qualities",
                "ymin": 0,
                "ylabel": "Percentage at base quality",
                "xlabel": "base quality",
                "colors": color_dict,
            },
        ],
        "id": "per_run_bq_hist",
        "title": "bases2fastq: Quality Histograms",
        "ylab": "Percentage",
    }
    plot_html = linegraph.plot(plot_content, pconfig=config)
    plot_name = "Per Run Base Quality Histogram"
    anchor = "bq_hist"
    description = "Histogram of base qualities"
    helptext = """
    This sections plot the base qualities histograms. 
    
    Per base quality histogram plots the distribution of each inidividual base qualities.
    
    Per read quality histogram plots the distribution of average base quality of each read
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content


def plot_base_quality_by_cycle(run_data, color_dict):
    # Prepare plot data for median BQ of each cycle

    r1r2_split = 0
    for s_name in run_data.keys():
        cycle_dict = dict()
        R1CycleNum = len(run_data[s_name]["Reads"][0]["Cycles"])
        r1r2_split = max(r1r2_split, R1CycleNum)

    median_dict = {}
    for s_name in run_data.keys():
        cycle_dict = dict()
        R1CycleNum = len(run_data[s_name]["Reads"][0]["Cycles"])
        R2CycleNum = len(run_data[s_name]["Reads"][0]["Cycles"])
        for cycle in run_data[s_name]["Reads"][0]["Cycles"]:
            cycle_no = cycle["Cycle"]
            cycle_dict.update({cycle_no: cycle["QualityScore50thPercentile"]})
        for cycle in run_data[s_name]["Reads"][1]["Cycles"]:
            cycle_no = int(cycle["Cycle"]) + r1r2_split
            cycle_dict.update({cycle_no: cycle["QualityScore50thPercentile"]})
        median_dict.update({s_name: cycle_dict})

    # Prepare plot data for mean BQ of each cycle
    mean_dict = {}
    for s_name in run_data.keys():
        ###Update each sample cycle info
        cycle_dict = dict()
        R1CycleNum = len(run_data[s_name]["Reads"][0]["Cycles"])
        R2CycleNum = len(run_data[s_name]["Reads"][0]["Cycles"])
        for cycle in run_data[s_name]["Reads"][0]["Cycles"]:
            cycle_no = cycle["Cycle"]
            cycle_dict.update({cycle_no: cycle["QualityScoreMean"]})
        for cycle in run_data[s_name]["Reads"][1]["Cycles"]:
            cycle_no = int(cycle["Cycle"]) + r1r2_split
            cycle_dict.update({cycle_no: cycle["QualityScoreMean"]})
        mean_dict.update({s_name: cycle_dict})

    # Prepare plot data for %Q30 of each cycle
    Q30_dict = {}
    for s_name in run_data.keys():
        ###Update each sample cycle info
        cycle_dict = dict()
        R1CycleNum = len(run_data[s_name]["Reads"][0]["Cycles"])
        R2CycleNum = len(run_data[s_name]["Reads"][0]["Cycles"])
        for cycle in run_data[s_name]["Reads"][0]["Cycles"]:
            cycle_no = cycle["Cycle"]
            cycle_dict.update({cycle_no: cycle["PercentQ30"]})
        for cycle in run_data[s_name]["Reads"][1]["Cycles"]:
            cycle_no = int(cycle["Cycle"]) + r1r2_split
            cycle_dict.update({cycle_no: cycle["PercentQ30"]})
        Q30_dict.update({s_name: cycle_dict})

    # Prepare plot data for %Q40 of each cycle
    Q40_dict = {}
    for s_name in run_data.keys():
        cycle_dict = dict()
        R1CycleNum = len(run_data[s_name]["Reads"][0]["Cycles"])
        R2CycleNum = len(run_data[s_name]["Reads"][0]["Cycles"])
        for cycle in run_data[s_name]["Reads"][0]["Cycles"]:
            cycle_no = cycle["Cycle"]
            cycle_dict.update({cycle_no: cycle["PercentQ40"]})
        for cycle in run_data[s_name]["Reads"][1]["Cycles"]:
            cycle_no = int(cycle["Cycle"]) + r1r2_split
            cycle_dict.update({cycle_no: cycle["PercentQ40"]})
        Q40_dict.update({s_name: cycle_dict})

    # aggregate plot data
    plot_content = [median_dict, mean_dict, Q30_dict, Q40_dict]
    config = {
        "data_labels": [
            {"name": "Median", "xlab": "cycle", "ylab": "Base Quality", "ymax": 50},
            {"name": "Mean", "ylab": "Quality", "ymax": 50},
            {"name": "%Q30", "xlab": "cycle", "ylab": "Percentage", "ymax": 100},
            {"name": "%Q40", "xlab": "cycle", "ylab": "Percentage", "ymax": 100},
        ],
        "xPlotLines": [{"color": "#FF0000", "width": 2, "value": r1r2_split, "dashStyle": "Dash"}],
        "colors": color_dict,
        "ymin": 0,
        "id": "per_run_quality_by_cycle",
        "title": "bases2fastq: Quality by cycles",
        "ylab": "QC",
    }
    plot_html = linegraph.plot(plot_content, pconfig=config)
    plot_name = "Quality statistics by cycle"
    anchor = "per_cycle_quality"
    description = "Per Run Base qualities by cycle"
    helptext = """
    This section plots the base qualities by each sequencing cycle. You can choose to show
    either median, mean, percentage Q30 or percentage Q40 of each cycle. Read 1 and Read 2 
    are separated by a red dashed line, such that we can visualize R1 and R2 qualities 
    in one plot.
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content
