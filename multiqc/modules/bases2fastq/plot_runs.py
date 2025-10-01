import math

from multiqc.plots import bargraph, linegraph, table
from multiqc import config
from natsort import natsorted
import random
import string

"""
Functions for plotting per run information of bases2fastq
"""


def generate_random_string():
    return "".join(random.choices(string.ascii_letters + string.digits, k=4))


def plot_run_stats(run_data, color_dict):
    """
    Plot a bar graph for polony numbers, Q30/Q40, index assignment rate and yields for each run
    """
    random_id = generate_random_string()
    run_names = list(run_data.keys())
    run_names.sort()
    num_polonies = dict()
    yields = dict()
    for run in run_names:
        # Index Assignment Polonies and Yields ###
        percent_assigned = run_data[run].get("PercentAssignedReads", 100.0)
        # percent_assigned = run_data[run]["PercentAssignedReads"]

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
        total_yield = run_data[run].get("TotalYield", run_data[run].get("AssignedYield", 300.0))
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
            {"name": "Yield (Gb)", "ylab": "Yield"},
        ],
        "cpswitch": True,
        "stacking": "normal",
        "id": f"run_metrics_bar_{random_id}",
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
    anchor = f"run_metrics_bar_{random_id}"
    description = "Bar plots of sequencing run yields. Please see individual run reports for details"
    helptext = """
    This section shows and compare the yield and index assignment rate of each sequencing run.\n\n
        - Number of Polonies: The total number of polonies that are calculated for the run.\n
        - Yield: The total yield of all assigned reads in gigabases.
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content


def _calculate_reads_eliminated(run_data) -> int:
    """
    Calculate the total number of reads eliminated during trimming.

    This function iterates over the lanes in the given run data and sums the
    difference between the number of polonies before trimming and after trimming.
    If required fields are missing, they are skipped.

    Args:
        run_data (dict): Dictionary containing sequencing run data with lane information.

    Returns:
        int: The total number of reads eliminated across all lanes.
    """
    reads_eliminated = 0
    if "Lanes" not in run_data:
        return reads_eliminated
    for lane in run_data["Lanes"]:
        if "NumPolonies" not in lane or "NumPoloniesBeforeTrimming" not in lane:
            continue
        reads_eliminated += lane["NumPoloniesBeforeTrimming"] - lane["NumPolonies"]

    return reads_eliminated


def tabulate_project_stats(run_data, color_dict):
    """
    Tabulate general information and statistics of each run
    """
    plot_content = dict()
    for s_name in run_data.keys():
        project = run_data[s_name]["Project"]
        run_project_name = f"{s_name} | {project}"
        run_stats = dict()
        run_stats.update({"num_polonies_run": int(run_data[s_name]["NumPolonies"])})
        run_stats.update({"yield_run": run_data[s_name]["AssignedYield"]})
        run_stats.update({"mean_base_quality_run": run_data[s_name]["QualityScoreMean"]})
        run_stats.update({"percent_q30_run": run_data[s_name]["PercentQ30"]})
        run_stats.update({"percent_q40_run": run_data[s_name]["PercentQ40"]})
        run_stats.update({"reads_eliminated": _calculate_reads_eliminated(run_data[s_name])})
        plot_content.update({run_project_name: run_stats})

    headers = {}
    headers["num_polonies_run"] = {
        "title": "# Polonies",
        "description": "The total number of polonies that are calculated for the run.",
        "min": 0,
        "scale": "RdYlGn",
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
        "title": "Assigned Yield (Gb)",
        "description": "The run yield based on assigned reads in gigabases",
        "scale": "Blues",
    }
    headers["mean_base_quality_run"] = {
        "title": "Quality Score Mean",
        "description": "Average base quality across Read 1 and Read 2",
        "min": 0,
        "scale": "Spectral",
    }
    headers["percent_q30_run"] = {
        "title": "Percent Q30",
        "description": "The percentage of ≥ Q30 Q scores for the project. This includes assigned and unassigned reads and excludes filtered reads and no calls.",
        "max": 100,
        "min": 0,
        "scale": "RdYlGn",
        "suffix": "%",
    }
    headers["percent_q40_run"] = {
        "title": "Percent Q40",
        "description": "The percentage of ≥ Q40 Q scores for the project. This includes assigned and unassigned reads and excludes filtered reads and no calls.",
        "max": 100,
        "min": 0,
        "scale": "RdYlGn",
        "suffix": "%",
    }
    headers["reads_eliminated"] = {
        "title": "Reads Eliminated",
        "description": "Number of reads eliminated.",
    }

    pconfig = {
        "title": "bases2fastq: General Sequencing (Project) QC metrics",
        "col1_header": "Run Name",
        "id": "project_run_metrics_table",
        "ylab": "QC",
    }

    project_header = ""
    run_keys = list(run_data.keys())
    if len(run_keys) > 1:
        project_header = "(Project) "
    elif len(run_keys) == 1:
        first_key = run_keys[0]
        project_header = f"{run_data[first_key]['Project']} | "
    plot_name = f"{project_header}Sequencing QC Metrics Table"
    plot_html = table.plot(plot_content, headers, pconfig=pconfig)
    anchor = "project_run_qc_metrics_table"
    description = "QC metrics per run, per project"
    helptext = """
    This section displays metrics that indicate the quality of each sequencing run: \n
        - Run Name: Unique identifier composed of (RunName)__(UUID), where (RunName) maps to the AVITI run name and (UUID) maps to the unique Bases2Fastq analysis result.\n
        - Number of Polonies: The total number of polonies that are calculated for the run.\n
        - Percentage Assigned Reads: The percentage of reads that are assigned to a sample.\n
        - Assigned Yield (Gb): The run yield that is based on assigned reads in gigabases.\n
        - Quality Score Mean: The mean Q score of base calls for the samples. This excludes filtered reads and no calls.\n
        - Percent Q30: The percentage of ≥ Q30 Q scores for the run. This includes assigned and unassigned reads and excludes filtered reads and no calls.\n
        - Percent Q40: The percentage of ≥ Q40 Q scores for the run. This includes assigned and unassigned reads and excludes filtered reads and no calls.\n
        - Reads Eliminated: Number of reads eliminated across lanes.\n
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
        run_stats.update({"percent_unexpected_index_pairs": run_data[s_name].get("PercentUnexpectedIndexPairs", 0.0)})
        run_stats.update({"yield_run": run_data[s_name]["AssignedYield"]})
        run_stats.update({"mean_base_quality_run": run_data[s_name]["QualityScoreMean"]})
        run_stats.update({"percent_q30_run": run_data[s_name]["PercentQ30"]})
        run_stats.update({"percent_q40_run": run_data[s_name]["PercentQ40"]})
        run_stats.update({"reads_eliminated": _calculate_reads_eliminated(run_data[s_name])})
        plot_content.update({s_name: run_stats})

    headers = {}
    headers["num_polonies_run"] = {
        "title": "# Polonies",
        "description": "The total number of polonies that are calculated for the run.)",
        "min": 0,
        "scale": "RdYlGn",
    }
    headers["percent_assigned_run"] = {
        "title": "% Assigned Reads",
        "description": "The percentage of reads assigned to sample(s)",
        "max": 100,
        "min": 0,
        "scale": "BuPu",
        "suffix": "%",
    }
    headers["percent_unexpected_index_pairs"] = {
        "title": "% Unexpected Index Pairs",
        "description": "The percentage of unexpected index pairs",
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
    headers["reads_eliminated"] = {
        "title": "Reads Eliminated",
        "description": "Number of reads eliminated.",
    }

    pconfig = {
        "title": "Bases2Fastq: General Sequencing Run QC metrics",
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
        - Reads Eliminated: Number of reads eliminated across lanes.\n
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content


def tabulate_manifest_stats(run_data, color_dict):
    """
    Tabulate general information and statistics of each run
    """
    random_id = generate_random_string()
    plot_content = dict()
    for s_name in run_data.keys():
        run_stats = dict()
        run_stats.update({"indexing": run_data[s_name]["Indexing"]})
        run_stats.update({"adapter_trim_type": run_data[s_name]["AdapterTrimType"]})
        run_stats.update({"min_read_length_r1": run_data[s_name]["R1AdapterMinimumTrimmedLength"]})
        run_stats.update({"min_read_length_r2": run_data[s_name]["R2AdapterMinimumTrimmedLength"]})
        plot_content.update({s_name: run_stats})

    headers = {}
    headers["indexing"] = {
        "title": "Indexing",
        "description": "Indexing scheme.",
        "scale": "RdYlGn",
    }
    headers["adapter_trim_type"] = {
        "title": "Adapter Trim Type",
        "description": "Adapter trimming method.",
    }
    headers["min_read_length_r1"] = {
        "title": "Minimum Read Length R1",
        "description": "Minimum read length for read R1.",
        "scale": "RdYlGn",
    }
    headers["min_read_length_r2"] = {
        "title": "Minimum Read Length R2",
        "description": "Minimum read length for read R1 (if applicable).",
        "scale": "RdYlGn",
    }

    pconfig = {
        "title": "Bases2Fastq: Run Manifest Metrics",
        "col1_header": "Run Name | Lane",
        # "id": f"run_manifest_metrics_table_{random_id}",
    }

    plot_name = "Run Manifest Table"
    plot_html = table.plot(plot_content, headers, pconfig=pconfig)
    anchor = f"run_manifest_metrics_table"
    description = "Run parameters used."
    helptext = """
    This section displays metrics that indicate the parameters used in the run: \n
        - Run Name | Lane: Unique identifier composed of (RunName)__(UUID) | (Lane), where (RunName) maps to the AVITI run name and (UUID) maps to the unique Bases2Fastq analysis result.\n
        - Indexing: Describes the indexing scheme.\n
        - Adapter Trim Type: Adapter trimming method.\n
        - Minimum Read Length R1/R2: Minumum read length after adapter trimming.\n
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content


def tabulate_index_assignment_stats(run_data, color_dict):
    """
    Tabulate general information and statistics of each run
    """
    plot_content = dict()
    sorted_run_data = natsorted(run_data.items(), key=lambda x: x[1]["SampleID"])
    for index, sample_data in enumerate(sorted_run_data, start=1):
        sample_data = sample_data[1]
        sample_index_stats = dict()
        sample_index_stats.update({"sample_name": sample_data["SampleID"]})
        sample_index_stats.update({"index_1": sample_data["Index1"]})
        sample_index_stats.update({"index_2": sample_data["Index2"]})
        sample_index_stats.update({"polonies": sample_data["SamplePolonyCounts"]})
        sample_index_stats.update({"polony_percentage": sample_data["PercentOfPolonies"]})
        plot_content.update({index: sample_index_stats})

    headers = {}
    headers["sample_name"] = {
        "title": "Sample Name",
        "description": "Sample Name (RunID + Sample ID).",
    }
    headers["index_1"] = {
        "title": "Index 1",
        "description": "Sample Index 1 (I1).",
    }
    headers["index_2"] = {
        "title": "Index 2",
        "description": "Sample Index 2 (I2).",
    }
    headers["polonies"] = {
        "title": "Polonies",
        "description": "Number of polonies assigned to sample.",
        "scale": "RdYlGn",
    }
    headers["polony_percentage"] = {
        "title": "Polony %",
        "description": "Percentage of total polonies assigned to this index combination.",
        "max": 100,
        "min": 0,
        "scale": "RdYlGn",
        "suffix": "%",
    }

    pconfig = {
        "title": "Bases2Fastq: Index Assignment Metrics",
        "col1_header": "Sample #",
        "id": "index_assignment_metrics",
    }

    plot_name = "Index Assignment Metrics"
    plot_html = table.plot(plot_content, headers, pconfig=pconfig)
    anchor = "index_assignment_metrics"
    description = "Index assignment metrics."
    helptext = """
    This section displays index assignment metrics including: \n
        - Sample Name: Sample identifier combining RunID and SampleID.\n
        - Index 1: Sample I1.\n
        - Index 2: Sample I2.\n
        - Polonies: Number of polonies assigned each sample.\n
        - Polony %: Percentage of total run's polonies assigned to each sample.\n
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content


def tabulate_unassigned_index_stats(run_data, color_dict):
    """
    Tabulate unassigned index metrics.

    run_data: Dictionary with unassigned index data including:
        - RunName
        - Lane
        - I1
        - I2
        - Polonies
        - % Polonies
    """

    headers = {}
    headers["Run Name"] = {
        "title": "Run Name",
        "description": "Run Name (Run ID + Analysis ID).",
    }
    headers["Lane"] = {
        "title": "Lane",
        "description": "Index Lane.",
    }
    headers["I1"] = {
        "title": "I1",
        "description": "Index 1.",
    }
    headers["I2"] = {
        "title": "I2",
        "description": "Index 2.",
    }
    headers["Polonies"] = {
        "title": "Polonies",
        "description": "Number of polonies assigned to indices.",
        "scale": "GnYlRd",
    }
    headers["% Polonies"] = {
        "title": "% Polonies",
        "description": "Percentage of total polonies assigned to this index combination.",
        "max": 100,
        "min": 0,
        "scale": "GnYlRd",
        "suffix": "%",
    }

    pconfig = {
        "title": "Bases2Fastq: Unassigned Indices Metrics",
        "col1_header": "Index #",
        "id": "index_unassignment_metrics",
    }

    plot_name = "Unassigned Indices Metrics"
    plot_html = table.plot(run_data, headers, pconfig=pconfig)
    anchor = "index_unassignment_metrics"
    description = "Index unassignment metrics."
    helptext = """
    This section displays index assignment metrics including: \n
        - Run Name: Run identifier. Built from Run ID and Analysis ID.\n
        - Lane: Lane number.\n
        - Index 1: Sample I1.\n
        - Index 2: Sample I2.\n
        - Polonies: Number of polonies assigned each index combination.\n
        - Polony %: Percentage of total run's polonies assigned to each index combination.\n
    """
    return plot_html, plot_name, anchor, description, helptext, run_data


def plot_lane_cycle_stats(run_data, color_dict):
    """
    Plot number of cycles per read and lane
    """
    plot_content = dict()
    for s_name in run_data.keys():
        if "Lanes" not in run_data[s_name]:
            continue
        for lane in run_data[s_name]["Lanes"]:
            if "Lane" not in lane or "Reads" not in lane:
                continue
            lane_stats = dict()
            lane_name = f"L{lane['Lane']}"
            run_name = f"{s_name} | {lane_name}"
            lane_stats[run_name] = {}
            for read in lane["Reads"]:
                if "Cycles" not in read or "Read" not in read:
                    continue
                read_name = read["Read"]
                num_cycles = len(read["Cycles"])
                lane_stats[run_name][read_name] = num_cycles
            plot_content.update(lane_stats)

    pconfig = {
        "title": "Bases2Fastq: Cycles Per Read Per Lane",
        "id": "project_cycles_per_read_per_lane",
        "ylab": "Read Cycles",
        "cpswitch": False,
        "subtitle": None,
    }

    plot_name = "Cycles Per Read Per Lane"
    plot_html = bargraph.plot(plot_content, pconfig=pconfig)
    anchor = "cycles_per_read_per_lane"
    description = "Number of sequencing cycles per read in each lane."
    helptext = """
    Shows the number of cycles used for each read in every flowcell lane. 
    Useful for confirming that read lengths match the expected sequencing setup across all lanes.
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
                "xlab": "Q Score",
                "colors": color_dict,
            },
            {
                "name": "Quality Per Read",
                "description": "Histogram of average read base quality",
                "ymin": 0,
                "ylabel": "Percentage of read quality",
                "xlab": "Q Score",
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
            {"name": "Mean Quality", "xlab": "cycle", "ylab": "Quality"},
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
    description = "Per run base qualities by cycle. Read 1 and Read 2 are separated by a red dashed line."
    helptext = """
    This section plots the base qualities by each instrument cycle.\n
    Choose between Median Quality, Mean Quality, Percent Q30 or Percentage Q40 per cycle.\n
    Read 1 and Read 2 are separated by a red dashed line.
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content
