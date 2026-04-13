from multiqc.plots import table

from .queries import (
    get_batch_counts,
    get_batch_density,
    get_cell_count,
    get_median_cell_diameter,
    get_percent_assigned,
    get_percent_confluency,
    get_percent_mismatch,
    get_percent_nucleated_cells,
    get_total_counts,
    get_total_density,
    get_percent_assigned_target_polony,
    get_percent_mismatch_target_polony,
    get_target_cell_metric_by_key,
    get_batch_extracellularratio,
)


def merge_well_dictionaries(dict_list):
    """
    Merge metric dictionaries to build a table of well metrics
    Input is list of dictionaries of the form
    {
        "{run_name} {well_location}" :
        {
            "{metric_name}" : "{metric_value}",
            ...
        },
        ...
    }
    Output will be a merged dictionary of the form
    {
        "{run_name} {well_location}" : {
            "{metric_name}" : "{metric_value}",
            ...
        },
        ...
    }
    """
    result = {}
    for data_dict in dict_list:
        for key in data_dict:
            key_dict = result.setdefault(key, {})
            for sub_key, sub_value in data_dict[key].items():
                key_dict[sub_key] = sub_value
    return result


def tabulate_wells(c2s_run_data):
    """
    Generate a table of well metrics from the cells2stats report
    """
    plot_content = merge_well_dictionaries(
        [
            get_cell_count(c2s_run_data),
            get_percent_confluency(c2s_run_data),
            get_percent_nucleated_cells(c2s_run_data),
            get_median_cell_diameter(c2s_run_data),
            get_total_density(c2s_run_data),
            get_total_counts(c2s_run_data),
        ]
    )
    headers = {}
    headers["cell_count"] = {
        "title": "# Cells",
        "description": "The number of cells in the well",
        "min": 0,
        "format": "{:d}",
        "scale": "GnBu",
    }

    headers["percent_confluency"] = {
        "title": "% Confluency",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["percent_nucleated_cells"] = {
        "title": "% Nucleated Cells",
        "description": "The percent of cells with a segmented nucleus",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["median_cell_diameter"] = {
        "title": "Median Cell Diameter",
        "description": "Median cell diameter for cells in the well",
        "min": 0,
        "scale": "GnBu",
        "suffix": "um",
        "format": "{:.1f}",
    }
    headers["total_density"] = {
        "title": "Assigned Counts K / mm2",
        "description": "Total density of assigned counts per mm2 of cell area across all barcoding batches",
        "min": 0,
        "scale": "GnBu",
    }
    headers["total_count"] = {
        "title": "Assigned Counts / Cell",
        "description": "Total average counts per cell across all barcoding batches",
        "min": 0,
        "scale": "GnBu",
    }

    pconfig = {
        "title": "cells2stats: Well cell paint and barcoding QC metrics",
        "col1_header": "Run / Well",
        "id": "well_metrics_table",
        "ylab": "QC",
    }

    plot_name = "Well Cell Paint and Barcoding QC Metrics Table"
    plot_html = table.plot(plot_content, headers, pconfig=pconfig)
    anchor = "well_metrics"
    description = "Table of cell paint and barocoding well QC metrics"
    helptext = """Provides cell paint and barocding metrics summarizing the performance of each well"""
    return plot_html, plot_name, anchor, description, helptext, plot_content


def tabulate_target_wells(c2s_run_data, target_site_name):
    """
    Generate a table of target well metrics from the cells2stats report
    """
    plot_content = merge_well_dictionaries(
        [
            get_percent_assigned_target_polony(c2s_run_data, target_site_name),
            get_percent_mismatch_target_polony(c2s_run_data, target_site_name),
            get_target_cell_metric_by_key(c2s_run_data, target_site_name, "PercentAssignedPureCells"),
            get_target_cell_metric_by_key(c2s_run_data, target_site_name, "PercentAssignedMixedCells"),
            get_target_cell_metric_by_key(c2s_run_data, target_site_name, "PercentUnassignedMixedCells"),
            get_target_cell_metric_by_key(c2s_run_data, target_site_name, "PercentUnassignedLowCountCells"),
            get_target_cell_metric_by_key(c2s_run_data, target_site_name, "AssignedCountsPerMM2", 1000),
            get_target_cell_metric_by_key(c2s_run_data, target_site_name, "MeanAssignedCountPerCell"),
            get_target_cell_metric_by_key(c2s_run_data, target_site_name, "MedianAbundantTargetCount"),
            get_target_cell_metric_by_key(c2s_run_data, target_site_name, "MeanUniqueTargetsPerCell"),
            get_target_cell_metric_by_key(c2s_run_data, target_site_name, "ExtraCellularRatio"),
            get_target_cell_metric_by_key(c2s_run_data, target_site_name, "PercentTargetDropout"),
        ]
    )
    headers = {}
    headers["percent_assigned"] = {
        "title": "% Assigned",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["percent_mismatch"] = {
        "title": "% Mismatch",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["PercentAssignedPureCells"] = {
        "title": "% Assigned Pure Cells",
        "description": "Percent of cells assigned to a target with one target present",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["PercentAssignedMixedCells"] = {
        "title": "% Assigned Mixed Cells",
        "description": "Percent of cells assigned to a target with greater than one target present",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["PercentUnassignedMixedCells"] = {
        "title": "% Unassigned Mixed Cells",
        "description": "Percent of cells not assigned to a target with greater than one target present",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["PercentUnassignedLowCountCells"] = {
        "title": "% Unassigned Low Count Cells",
        "description": "Percent of cells not assigned to a target with zero or one targets present",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }
    headers["AssignedCountsPerMM2"] = {
        "title": "Assigned Counts K / mm2",
        "description": "Total density of assigned target counts per mm2 of cell area",
        "min": 0,
        "scale": "GnBu",
        "suffix": "K",
    }
    headers["MeanAssignedCountPerCell"] = {
        "title": "Assigned Counts / Cell",
        "description": "Average assigned target counts per cell",
        "min": 0,
        "scale": "GnBu",
        "suffix": "",
    }
    headers["MedianAbundantTargetCount"] = {
        "title": "Median Abundant Target Count",
        "description": "Median most abundant target count for cells in the well",
        "min": 0,
        "scale": "GnBu",
        "suffix": "",
    }
    headers["MeanUniqueTargetsPerCell"] = {
        "title": "Mean Unique Targets / Cell",
        "description": "Mean number of unique targets per cell in the well",
        "min": 0,
        "scale": "GnBu",
        "suffix": "",
    }
    headers["ExtraCellularRatio"] = {
        "title": "Extracellular Ratio",
        "description": "Ratio of density of extra-cellular counts to density of intra-cellular counts",
        "min": 0,
        "scale": "GnBu",
        "suffix": "",
        "format": "{:.2f}",
    }
    headers["PercentTargetDropout"] = {
        "title": "% Target Dropout",
        "description": "Percent of targets not assigned to any cell",
        "scale": "GnBu",
        "max": 100,
        "min": 0,
        "suffix": "%",
    }

    pconfig = {
        "title": f"cells2stats: {target_site_name} target site well QC metrics",
        "col1_header": "Target Group / Well",
        "id": f"{target_site_name}_target_metrics_table",
        "ylab": "QC",
    }

    plot_name = f"{target_site_name} Target Site Well QC Metrics Table"
    plot_html = table.plot(plot_content, headers, pconfig=pconfig)
    anchor = f"{target_site_name}_target_well_metrics"
    description = f"Table of well QC metrics for target site {target_site_name}"
    helptext = f"Provides metrics summarizing the performance of each well for target site {target_site_name}"
    return plot_html, plot_name, anchor, description, helptext, plot_content


def merge_batch_dictionaries(dict_list, metric_names):
    """
    Merge input dictionaries to build a table of batch metrics
    Input is list of dictionaries of the form
    {
        "{run_name} {well_location}" :
        {
            "{batch_name}" : "{batch_value}",
            ...
        },
        ...
    }

    Output will be a merged dictionary of the form
    {
        "{run name} {well_location} {batch_name}" : {
            "{metric_name}: "{metric_value}",
            ...
        },
        ...
    }
    """
    assert len(dict_list) == len(metric_names), "Input dictionary list and metric names must be the same length"
    result = {}
    for data_dict, metric_name in zip(dict_list, metric_names):
        for well_name in data_dict:
            well_dict = data_dict[well_name]
            for batch_name in well_dict:
                result.setdefault(f"{well_name} {batch_name}", {})[metric_name] = well_dict[batch_name]
    return result


def tabulate_batches(c2s_run_data):
    """
    Generate a table of barcoding batch metrics from the cells2stats report
    """
    plot_content = merge_batch_dictionaries(
        [
            get_batch_density(c2s_run_data),
            get_batch_counts(c2s_run_data),
            get_percent_assigned(c2s_run_data),
            get_percent_mismatch(c2s_run_data),
            get_batch_extracellularratio(c2s_run_data),
        ],
        ["batch_density", "batch_count", "percent_assigned", "percent_mismatch", "extracellular_ratio"],
    )

    headers = {}
    headers["batch_density"] = {
        "title": "Assigned Counts K / mm2",
        "description": "Assigned counts per mm2 for each barcoding batch in the well",
        "min": 0,
        "scale": "GnBu",
    }
    headers["batch_count"] = {
        "title": "Assigned Counts / Cell",
        "description": "Average assigned counts per cell for each barcoding batch in the well",
        "scale": "GnBu",
        "min": 0,
    }

    headers["percent_assigned"] = {
        "title": "% Assigned",
        "description": "Percent of polonies assigned to a barcode",
        "scale": "GnBu",
        "min": 0,
    }
    headers["percent_mismatch"] = {
        "title": "% Mismatch",
        "description": "Percent of assigned polonies assigned with a mismatch",
        "scale": "GnBu",
        "min": 0,
    }
    headers["extracellular_ratio"] = {
        "title": "Extracellular Ratio",
        "description": "Ratio of density of extra-cellular counts to density of intra-cellular counts",
        "min": 0,
        "scale": "GnBu",
        "format": "{:.2f}",
    }

    pconfig = {
        "title": "cells2stats: Barcodin batch QC metrics",
        "col1_header": "Run / Well / Batch",
        "id": "batch_metrics_table",
        "ylab": "QC",
    }

    plot_name = "Barcoding Batch QC Metrics Table"
    plot_html = table.plot(plot_content, headers, pconfig=pconfig)
    anchor = "batch_metrics"
    description = "Table of barcoding batch QC metrics"
    helptext = """Provides overall metrics summarizing the performance of each barcoding batch per well"""
    return plot_html, plot_name, anchor, description, helptext, plot_content
