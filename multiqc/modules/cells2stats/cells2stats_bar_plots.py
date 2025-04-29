from multiqc.utils import mqc_colour
from multiqc.plots import bargraph

from .utils import summarize_batch_names
from .queries import *

def plot_barcoding(c2s_run_data):
    """"
    Generate plots related to barcoding performance metrics from the cells2stats report
    """

    batch_names = summarize_batch_names(c2s_run_data)    
    plot_content = [get_percent_assigned(c2s_run_data), get_percent_mismatch(c2s_run_data),]
    pconfig = {
        "data_labels": [
            {"name": "Assigned Reads", "ylab": "Percent of Assigned Reads"},
            {"name" : "Mismatches", "ylab" : "Percent Mismatch"}

        ],
        "cpswitch": False,
        "id": "barcoding_bar",
        "stacking": "group",  
        "title": "cells2stats: Barcoding QC metrics plot",
        "ylab": "QC"
    }

    scale = mqc_colour.mqc_colour_scale("GnBu", 0, len(batch_names))
    
    cat = {}
    for i, batch_name in enumerate(batch_names):
        cat[batch_name] = {"name": batch_name, "color": scale.get_colour(i, lighten=1)}
    
    cats = [
        cat,
        cat
    ]   

    plot_name = "Barcoding Metrics"
    plot_html = bargraph.plot(plot_content, cats, pconfig=pconfig)
    anchor = "well_barcoding_plot"
    description = "Bar plots of barcoding metrics"
    helptext = """PLot percent of assigned reads and percent of reads assigned with mismatch for each well in the run."""
    
    return plot_html, plot_name, anchor, description, helptext, plot_content


def plot_cell_assignment(c2s_run_data):
    """"
    Generate plots related to cell assignment performance metrics from the cells2stats report
    """
        
    batch_names = summarize_batch_names(c2s_run_data)
    plot_content = [get_total_density(c2s_run_data), get_batch_density(c2s_run_data), get_total_counts(c2s_run_data), get_batch_counts(c2s_run_data),]
    pconfig = {
        "data_labels": [
            {"name": "Total Density", "ylab": "Assigned Counts K / mm2"},
            {"name": "Batch Density", "ylab": "Assigned Counts K / mm2"},
            {"name": "Total Counts", "ylab": "Average Assigned Counts / Cell"},
            {"name": "Batch Counts", "ylab": "Average Assigned Counts / Cell"},
        ],
        "cpswitch": False,
        "id": "cell_assignment_bar",
        "stacking": "group",  
        "title": "cells2stats: Barcoding QC metrics plot",
        "ylab": "QC"
    }

    scale = mqc_colour.mqc_colour_scale("GnBu", 0, len(batch_names))
    
    cat = {}
    for i, batch_name in enumerate(batch_names):
        cat[batch_name] = {"name": batch_name, "color": scale.get_colour(i, lighten=1)}
    
    cats = [
        {"total_density" : {"name": "Total Density"} },
        cat,
        {"total_count" : {"name": "Total Counts"} },
        cat
    ]   

    plot_name = "Cell Assignment Metrics"
    plot_html = bargraph.plot(plot_content, cats, pconfig=pconfig)
    anchor = "well_assignment_plot"
    description = "Bar plots of cell assignment metrics"
    helptext = """Plot density and absolute of assigned counts per batch and across all batches."""
    
    return plot_html, plot_name, anchor, description, helptext, plot_content


def plot_cell_segmentation(c2s_run_data):
    """"
    Generate plots related to cell segmentation metrics from the cells2stats report
    """
    plot_content = []
    plot_content.append(get_cell_count(c2s_run_data))
    plot_content.append(get_percent_confluency(c2s_run_data))
    plot_content.append(get_percent_nucleated_cells(c2s_run_data))
    plot_content.append(get_median_cell_diameter(c2s_run_data))

    pconfig = {
        "data_labels": [
            {"name": "Cell Count", "ylab": "Number of Cells", "format": "{d}"},
            {"name": "Confluency", "ylab": "Percent Confluency"},
            {"name": "Nucleated Cells", "ylab": "Percent Nucleated Cells"},
            {"name": "Cell Diameter", "ylab": "Median Cell Diameter (um)"},

        ],
        "cpswitch": False,
        "id": "cell_segmentation_bar",
        "title": "cells2stats: Cell segmentation QC metrics plot",
        "ylab": "QC",
    }

    cats = [
        {
            "cell_count": {"name": "Cell Count"}
        },
        {
            "percent_confluency": {"name": "Confluency"}
        },
        {
            "percent_nucleated_cells": {"name": "Nucleated Cells"}
        },
        {
            "median_cell_diameter": {"name": "Median Cell Diameter (um)"}
        }
    ]   

    plot_name = "Cell Segmentation Metrics"
    plot_html = bargraph.plot(plot_content, cats, pconfig=pconfig)
    anchor = "well_segmentation_plot"
    description = "Bar plots of cell segmentation metrics"
    helptext = """Plot cell counts, confluency, nucleated cells, and median cell diameter for each well in the run."""
    
    return plot_html, plot_name, anchor, description, helptext, plot_content

def plot_controls(c2s_run_data):
    """"
    Generate plots related to control metrics from the cells2stats report
    """
    
    batch_names = summarize_batch_names(c2s_run_data)
    
    controls = set()
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells",[]):
            for batch_data in well_data.get("Batches", []):
                batch_name = batch_data.get("BatchName", "")
                if batch_name != "" and not batch_name.startswith("CP"):
                    for control_target in batch_data.get("ControlTargets", []):
                        control_name = control_target.get("TargetName", "")
                        if control_name != "":
                            controls.add(control_name)
    controls = sorted(controls)[:8]

    plot_content = []
    for control in controls:
        control_content = {}
        for run_name in c2s_run_data:
            run_data = c2s_run_data[run_name]
            for well_data in run_data.get("CytoStats", {}).get("Wells", []):
                well_location = well_data.get("WellLocation", "")
                if well_location != "":
                    for batch_name in batch_names:
                        batch_data = find_entry(well_data.get("Batches", []), "BatchName", batch_name, {})
                        control_data = find_entry(batch_data.get("ControlTargets", []), "TargetName", control, {})
                        val = json_decode_float(control_data.get("AssignedCountPerMM2", float("nan")))
                        if not is_nan(val):
                            control_content.setdefault(f"{run_name} {well_location}",{})[batch_name] = val
        plot_content.append(control_content)

    pconfig = {
        "data_labels": [
            {"name" : control, "ylab" : "Assigned Counts K / mm2"} for control in controls],
        "cpswitch": False,
        "id": "controls_bar",
        "stacking": "group",  
        "title": "cells2stats: Control targets plot",
        "ylab": "QC"
    }

    scale = mqc_colour.mqc_colour_scale("GnBu", 0, len(batch_names))
    
    cat = {}
    for i, batch_name in enumerate(batch_names):
        cat[batch_name] = {"name": batch_name, "color": scale.get_colour(i, lighten=1)}
    
    cats = [cat,] * len(plot_content)

    plot_name = "Control Targets"
    plot_html = bargraph.plot(plot_content, cats, pconfig=pconfig)
    anchor = "well_control_plot"
    description = "Bar plots of control targets"
    helptext = """Plot density of assigned counts per batch for controls."""
    
    return plot_html, plot_name, anchor, description, helptext, plot_content
