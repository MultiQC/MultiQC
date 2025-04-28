from .utils import json_decode_int, json_decode_float, is_nan, summarize_batch_names, find_entry

def get_cell_count(c2s_run_data):
    """
    Get the cell count for each well in the run"""
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            result.setdefault(f"{run_name} {well_data['WellLocation']}",{})["cell_count"] = json_decode_int(well_data.get("CellCount", 0))
    return result


def get_percent_confluency(c2s_run_data):
    """
    Get the percent confluency for each well in the run
    """
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            val = json_decode_float(well_data.get("PercentConfluency", float("nan")))
            if not is_nan(val):
                result.setdefault(f"{run_name} {well_data['WellLocation']}",{})["percent_confluency"] = val
    return result


def get_percent_nucleated_cells(c2s_run_data):
    """
    Get the percent nucleated cells for each well in the run
    """
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            val = json_decode_float(well_data.get("PercentNucleatedCells", float("nan")))
            if not is_nan(val):
                result.setdefault(f"{run_name} {well_data['WellLocation']}",{})["percent_nucleated_cells"] = val
    return result


def get_median_cell_diameter(c2s_run_data):
    """
    Get the median cell diameter for each well in the run
    """
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            val = json_decode_float(well_data.get("MedianCellDiameter", float("nan")))
            if not is_nan(val):
                result.setdefault(f"{run_name} {well_data['WellLocation']}",{})["median_cell_diameter"] = val
    return result


def get_total_density(c2s_run_data):
    """
    Get the total density for each well in the run
    """
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            well_location = well_data.get("WellLocation", "")
            if well_location != "":
                val = json_decode_float(well_data.get("AssignedCountsPerMM2", float("nan")))
                if not is_nan(val):
                    result.setdefault(f"{run_name} {well_location}",{})["total_density"] = val / 1000.0
    return result


def get_total_counts(c2s_run_data):
    """
    Get the average total counts per cell for each well in the run
    """
    batch_names = summarize_batch_names(c2s_run_data)
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            well_location = well_data.get("WellLocation", "")
            if well_location != "":
                total_count = 0
                for batch_name in batch_names:
                    batch_data = find_entry(well_data.get("Batches", []), "BatchName", batch_name, {})
                    val = json_decode_float(batch_data.get("MeanAssignedCountPerCell", float("nan")))
                    if not is_nan(val):
                        total_count += val
                result.setdefault(f"{run_name} {well_location}",{})["total_count"] = total_count
    return result


def get_batch_density(c2s_run_data):
    """
    Get the average density for each batch in each well in the run"""
    batch_names = summarize_batch_names(c2s_run_data)
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            well_location = well_data.get("WellLocation", "")
            if well_location != "":
                for batch_name in batch_names:
                    batch_data = find_entry(well_data.get("Batches", []), "BatchName", batch_name, {})                  
                    val = json_decode_float(batch_data.get("AssignedCountsPerMM2", float("nan")))
                    if not is_nan(val):
                        result.setdefault(f"{run_name} {well_location}",{})[batch_name] = val / 1000.0
    return result


def get_batch_counts(c2s_run_data):
    """
    Get the average counts per cell for each batch in each well in the run
    """
    batch_names = summarize_batch_names(c2s_run_data)
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            well_location = well_data.get("WellLocation", "")
            if well_location != "":
                for batch_name in batch_names:
                    batch_data = find_entry(well_data.get("Batches", []), "BatchName", batch_name, {})                  
                    val = json_decode_float(batch_data.get("MeanAssignedCountPerCell", float("nan")))
                    if not is_nan(val):
                        result.setdefault(f"{run_name} {well_location}",{})[batch_name] = val
    return result


def get_percent_assigned(c2s_run_data):
    """
    Get the percent assigned reads for each batch in each well in the run
    """
    batch_names = summarize_batch_names(c2s_run_data)
    result = {}
    for run_name in c2s_run_data:
            run_data = c2s_run_data[run_name]
            for batch_name in batch_names:
                batch_entry = find_entry(run_data.get("DemuxStats", {}).get("Batches", []), "BatchName", batch_name, {})
                for well_data in batch_entry.get("Wells", []):
                    well_location = well_data.get("WellLocation", "")
                    if well_location != "":
                        val = json_decode_float(well_data.get("PercentAssignedReads", float("nan")))
                        if not is_nan(val):
                            result.setdefault(f"{run_name} {well_location}",{})[batch_name] = val
    return result


def get_percent_mismatch(c2s_run_data):
    """
    Get the percent mismatch for each batch in each well in the run
    """
    batch_names = summarize_batch_names(c2s_run_data)
    result = {}
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for batch_name in batch_names:
            batch_entry = find_entry(run_data.get("DemuxStats", {}).get("Batches", []), "BatchName", batch_name, {})
            for well_data in batch_entry.get("Wells", []):
                well_location = well_data.get("WellLocation", "")
                if well_location != "":
                    val = json_decode_float(well_data.get("PercentMismatch", float("nan")))
                    if not is_nan(val):
                        result.setdefault(f"{run_name} {well_location}",{})[batch_name] = val
    return result