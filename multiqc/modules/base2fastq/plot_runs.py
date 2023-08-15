from multiqc.plots import bargraph, table, linegraph, scatter
import numpy as np

"""
Functions for plotting per run information of base2fastq
"""


def plot_run_stats(runData, colorDict):
    """
    plot bargraph for polony numbers, Q30/Q40, index assignment rate and yields for each run
    """
    runNames = list(runData.keys())
    runNames.sort()
    NumPolonies = dict()
    PercentQs = dict()
    PercentAssigned = dict()
    Yields = dict()
    for run in runNames:
        NumPolonies.update({run: {"Number of Polonies": runData[run]["NumPolonies"]}})
        PercentQs.update({run: {"percent Q30": runData[run]["PercentQ30"], "percent Q40": runData[run]["PercentQ40"]}})
        PercentAssigned.update({run: {"Successful assigned reads": runData[run]["PercentAssignedReads"]}})
        Yields.update({run: {"Assigned yield": runData[run]["AssignedYield"]}})
    plotContent = [NumPolonies, PercentQs, PercentAssigned, Yields]
    config = {
        "data_labels": [
            {"name": "Polony Numbers", "ylab": "Number of polonies", "format": "{d}"},
            {"name": "Base Qualities", "ylab": "Percentage"},
            {"name": "Percentage assigned", "ylab": "Percentage"},
            {"name": "Data Yield", "ylab": "Gb"},
        ],
        "cpswitch": False,
        "stacking": None,
        "description": "Bar plots to compare some general statistics of all found sequencing runs.",
    }
    plotName = "Sequencing Run QC metrics"
    plotHtml = bargraph.plot(plotContent, pconfig=config)
    anchor = "run_qc_metrics_plot"
    description = "bar plots of general QC metrics"
    helptext = """
    This section shows and compare some metrics that indicate the quality of each sequencing run. 
    
    Polony numbers: total number of polonies. Each polony will yield one read1 and one read2
    Base Qualities: percentage of bases that has >=30 (blue) or >=40 (black) base qualities
    Percentage assinged: percentage of reads that has been assigned to any sample
    Data Yield: the total volume of data that has been assigned to any sample
    """
    return plotHtml, plotName, anchor, description, helptext


def tabulate_run_stats(runData, colorDict):
    """
    Tabulate general information and statistics of each run
    """
    plotContent = dict()
    for s_name in runData.keys():
        runStats = dict()
        runStats.update({"#Polonies": runData[s_name]["NumPolonies"]})
        runStats.update({"Percentage Assigned": runData[s_name]["PercentAssignedReads"]})
        runStats.update({"Percentage Q30": runData[s_name]["PercentQ30"]})
        runStats.update({"PercentageQ40": runData[s_name]["PercentQ40"]})
        runStats.update({"Average Base Quality": runData[s_name]["QualityScoreMean"]})
        runStats.update({"Yield(Gb)": runData[s_name]["AssignedYield"]})
        runStats.update({"Read1 Average Length": runData[s_name]["Reads"][0]["MeanReadLength"]})
        runStats.update({"Read2 Average Length": runData[s_name]["Reads"][1]["MeanReadLength"]})
        plotContent.update({s_name: runStats})
    headers = {header: {"title": header} for header in runStats.keys()}
    headers["#Polonies"].update({"format": "{d}"})
    config = {"descriptions": "Table of per sequencing run key informations", "col1_header": "Run Name"}
    plotName = "Sequencing Run basic informations"
    plotHtml = table.plot(plotContent, headers, pconfig=config)
    anchor = "run_qc_metrics_table"
    description = "table of general QC metrics"
    helptext = """
    This section shows numbers of some metrics that indicate the quality of each sequencing run. 
    \n
    Polony numbers: total number of polonies. Each polony yields one read1 and one read2
    \n
    Percentage assinged: percentage of reads that has been assigned to any sample
    \n
    Percentage Q30: percentage of bases that has >=30 base quality (error rate <= 10^-3) 
    \n
    Percentage Q40: percentage of bases that has >=40 base quality (error rate <= 10^-4)
    \n
    Average Base Quality: average base qualities across all bases
    \n
    Yield: total volume of data that has been assigned to any sample
    \n
    Read1 Average Length: the average length of Read 1
    \n
    Read2 Average Length: the average length of Read 2
    \n
    """
    return plotHtml, plotName, anchor, description, helptext


def plot_base_quality_hist(runData, colorDict):
    # Prepare plot data for per base BQ histogram
    bqHistDict = dict()
    for s_name in runData.keys():
        R1BaseQualityCounts = runData[s_name]["Reads"][0]["QualityScoreHistogram"]
        R2BaseQualityCounts = runData[s_name]["Reads"][1]["QualityScoreHistogram"]
        R1R2BaseQualityCounts = [r1 + r2 for r1, r2 in zip(R1BaseQualityCounts, R2BaseQualityCounts)]
        totalBases = sum(R1R2BaseQualityCounts)
        bqHistDict.update({s_name: {}})
        for quality in range(0, len(R1R2BaseQualityCounts)):
            bqHistDict[s_name].update({quality: R1R2BaseQualityCounts[quality] / totalBases * 100})

    # Prepare plot data for per read average BQ histogram
    perReadQualityHistDict = dict()
    for s_name in runData.keys():
        R1QualityCounts = runData[s_name]["Reads"][0]["PerReadMeanQualityScoreHistogram"]
        R2QualityCounts = runData[s_name]["Reads"][1]["PerReadMeanQualityScoreHistogram"]
        totalReads = runData[s_name]["NumPolonies"] * 2
        R1R2QualityCounts = [r1 + r2 for r1, r2 in zip(R1QualityCounts, R2QualityCounts)]
        perReadQualityHistDict.update({s_name: {}})
        for meanQuality in range(0, len(R1R2QualityCounts)):
            perReadQualityHistDict[s_name].update({meanQuality: R1R2QualityCounts[meanQuality] / totalReads * 100})
    plotContent = [bqHistDict, perReadQualityHistDict]

    # Config for switching dataset
    config = {
        "data_labels": [
            {
                "name": "per base quality histogram",
                "descriptions": "Histogram of bases qualities",
                "ymin": 0,
                "ylabel": "Percentage at base quality",
                "xlabel": "base quality",
                "colors": colorDict,
            },
            {
                "name": "per read quality histogram",
                "descriptions": "Histogram of read average base qualities",
                "ymin": 0,
                "ylabel": "Percentage at base quality",
                "xlabel": "base quality",
                "colors": colorDict,
            },
        ],
    }
    plotHtml = linegraph.plot(plotContent, pconfig=config)
    plotName = "Base Quality Histogram"
    anchor = "bq_hist"
    description = "Histogram of base qualities"
    helptext = """
    This sections plot the base qualities histograms. 
    
    Per base quality histogram plots the distribution of each inidividual base qualities.
    
    Per read quality histogram plots the distribution of average base quality of each read
    """
    return plotHtml, plotName, anchor, description, helptext


def plot_base_quality_by_cycle(runData, colorDict):
    # Prepare plot data for median BQ of each cycle
    medianDict = {}
    for s_name in runData.keys():
        cycleDict = dict()
        R1CycleNum = len(runData[s_name]["Reads"][0]["Cycles"])
        R2CycleNum = len(runData[s_name]["Reads"][0]["Cycles"])
        for cycle in runData[s_name]["Reads"][0]["Cycles"]:
            cycleNo = cycle["Cycle"]
            cycleDict.update({cycleNo: cycle["QualityScore50thPercentile"]})
        for cycle in runData[s_name]["Reads"][1]["Cycles"]:
            cycleNo = int(cycle["Cycle"]) + R1CycleNum
            cycleDict.update({cycleNo: cycle["QualityScore50thPercentile"]})
        medianDict.update({s_name: cycleDict})

    # Prepare plot data for mean BQ of each cycle
    meanDict = {}
    for s_name in runData.keys():
        ###Update each sample cycle info
        cycleDict = dict()
        R1CycleNum = len(runData[s_name]["Reads"][0]["Cycles"])
        R2CycleNum = len(runData[s_name]["Reads"][0]["Cycles"])
        for cycle in runData[s_name]["Reads"][0]["Cycles"]:
            cycleNo = cycle["Cycle"]
            cycleDict.update({cycleNo: cycle["QualityScoreMean"]})
        for cycle in runData[s_name]["Reads"][1]["Cycles"]:
            cycleNo = int(cycle["Cycle"]) + R1CycleNum
            cycleDict.update({cycleNo: cycle["QualityScoreMean"]})
        meanDict.update({s_name: cycleDict})

    # Prepare plot data for %Q30 of each cycle
    Q30Dict = {}
    for s_name in runData.keys():
        ###Update each sample cycle info
        cycleDict = dict()
        R1CycleNum = len(runData[s_name]["Reads"][0]["Cycles"])
        R2CycleNum = len(runData[s_name]["Reads"][0]["Cycles"])
        for cycle in runData[s_name]["Reads"][0]["Cycles"]:
            cycleNo = cycle["Cycle"]
            cycleDict.update({cycleNo: cycle["PercentQ30"]})
        for cycle in runData[s_name]["Reads"][1]["Cycles"]:
            cycleNo = int(cycle["Cycle"]) + R1CycleNum
            cycleDict.update({cycleNo: cycle["PercentQ30"]})
        Q30Dict.update({s_name: cycleDict})

    # Prepare plot data for %Q40 of each cycle
    Q40Dict = {}
    for s_name in runData.keys():
        cycleDict = dict()
        R1CycleNum = len(runData[s_name]["Reads"][0]["Cycles"])
        R2CycleNum = len(runData[s_name]["Reads"][0]["Cycles"])
        for cycle in runData[s_name]["Reads"][0]["Cycles"]:
            cycleNo = cycle["Cycle"]
            cycleDict.update({cycleNo: cycle["PercentQ40"]})
        for cycle in runData[s_name]["Reads"][1]["Cycles"]:
            cycleNo = int(cycle["Cycle"]) + R1CycleNum
            cycleDict.update({cycleNo: cycle["PercentQ40"]})
        Q40Dict.update({s_name: cycleDict})

    # aggregate plot data
    plotContent = [medianDict, meanDict, Q30Dict, Q40Dict]
    config = {
        "data_labels": [
            {"name": "Median", "xlab": "cycle", "ylab": "Base Quality", "ymax": 50},
            {"name": "Mean", "ylab": "Quality", "ymax": 50},
            {"name": "%Q30", "xlab": "cycle", "ylab": "Percentage", "ymax": 100},
            {"name": "%Q40", "xlab": "cycle", "ylab": "Percentage", "ymax": 100},
        ],
        "xPlotLines": [{"color": "#FF0000", "width": 2, "value": R1CycleNum, "dashStyle": "Dash"}],
        "colors": colorDict,
        "ymin": 0,
    }
    plotHtml = linegraph.plot(plotContent, pconfig=config)
    plotName = "Quality statistics by cycle"
    anchor = "per_cycle_quality"
    description = "Base qualities by cycle"
    helptext = """
    This section plots the base qualities by each sequencing cycle. You can choose to show
    either median, mean, percentage Q30 or percentage Q40 of each cycle. Read 1 and Read 2 
    are separated by a red dashed line, such that we can visualize R1 and R2 qualities 
    in one plot.
    """
    return plotHtml, plotName, anchor, description, helptext
