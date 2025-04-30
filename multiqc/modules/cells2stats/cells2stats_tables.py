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
        "scale": "GnBu",
        "format": "{d}",
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
        "description": "The percentage of cells with a segmented nucleus",
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
    }
    headers["total_density"] = {
        "title": "Assigned Counts K / mm2",
        "description": "Total density of assigned counts per mm2 of cell area across all batches",
        "min": 0,
        "scale": "GnBu",
    }
    headers["total_count"] = {
        "title": "Assigned Counts / Cell",
        "description": "Total average counts per cell across all batches",
        "min": 0,
        "scale": "GnBu",
    }

    pconfig = {
        "title": "cells2stats: Well QC metrics",
        "col1_header": "Run / Well",
        "id": "well_metrics_table",
        "ylab": "QC",
    }

    plot_name = "Well QC metrics table"
    plot_html = table.plot(plot_content, headers, pconfig=pconfig)
    anchor = "well_metrics_table"
    description = "Table of general well QC metrics"
    helptext = """Provides overall metrics summarizing the performance of each well"""
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
    Generate a table of batch metrics from the cells2stats report
    """
    plot_content = merge_batch_dictionaries(
        [
            get_batch_density(c2s_run_data),
            get_batch_counts(c2s_run_data),
            get_percent_assigned(c2s_run_data),
            get_percent_mismatch(c2s_run_data),
        ],
        ["batch_density", "batch_count", "percent_assigned", "percent_mismatch"],
    )

    headers = {}
    headers["batch_density"] = {
        "title": "Assigned Counts K / mm2",
        "description": "Assigned counts per mm2 for each batch in the well",
        "min": 0,
        "scale": "GnBu",
    }
    headers["batch_count"] = {
        "title": "Assigned Counts / Cell",
        "description": "Average assigned counts per cell for each batch in the well",
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

    pconfig = {
        "title": "cells2stats: Batch QC metrics",
        "col1_header": "Batch",
        "id": "batch_metrics_table",
        "ylab": "QC",
    }

    plot_name = "Batch QC metrics table"
    plot_html = table.plot(plot_content, headers, pconfig=pconfig)
    anchor = "batch_metrics_table"
    description = "Table of general batch QC metrics"
    helptext = """Provides overall metrics summarizing the performance of each batch per well"""
    return plot_html, plot_name, anchor, description, helptext, plot_content
