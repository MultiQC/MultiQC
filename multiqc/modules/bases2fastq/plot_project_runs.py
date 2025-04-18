from multiqc.plots import table
from multiqc import config

"""
Functions for plotting per run information of bases2fastq
"""


def tabulate_project_run_stats(run_data, color_dict):
    """
    Tabulate general information and statistics of each run
    """
    plot_content = dict()
    for s_name in run_data.keys():
        run_stats = dict()
        run_stats.update({"num_polonies_run": int(run_data[s_name]["NumPolonies"])})
        run_stats.update({"yield_run": run_data[s_name]["AssignedYield"]})
        run_stats.update({"mean_base_quality_run": run_data[s_name]["QualityScoreMean"]})
        run_stats.update({"percent_q30_run": run_data[s_name]["PercentQ30"]})
        run_stats.update({"percent_q40_run": run_data[s_name]["PercentQ40"]})
        plot_content.update({s_name: run_stats})

    headers = {}
    headers["num_polonies_run"] = {
        "title": f"# Polonies ({config.base_count_prefix})",
        "description": f"The total number of polonies that are calculated for the run ({config.base_count_desc})",
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

    pconfig = {
        "title": "bases2fastq: General Sequencing (Project) QC metrics",
        "col1_header": "Run Name",
        "id": "project_run_metrics_table",
        "ylab": "QC",
    }

    plot_name = "(Project) Sequencing QC metrics table"
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
    """
    return plot_html, plot_name, anchor, description, helptext, plot_content
