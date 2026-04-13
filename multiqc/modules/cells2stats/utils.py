def find_entry(json_list, key, value, default_value):
    """
    Find an entry in a list of dictionaries by key and value
    """
    for entry in json_list:
        if entry.get(key, "") == value:
            return entry
    return default_value


def json_decode_float(value):
    """
    Decode a JSON float value, converting it to NaN if it is -999
    """
    value = float(value)
    if abs(value + 999) < 0.00001:
        return float("nan")
    return value


def json_decode_int(value):
    """
    Decode a JSON int value, converting it to 0 if it is -999
    """
    value = int(value)
    if abs(value + 999) < 0.00001:
        return 0
    return value


def is_nan(value):
    """
    Check if a value is NaN
    """
    return value != value


def summarize_batch_names(c2s_run_data):
    """
    Generate a list of barcoding batch names from the cells2stats report
    """
    batch_names = set()
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for well_data in run_data.get("CytoStats", {}).get("Wells", []):
            for batch_data in well_data.get("Batches", []):
                batch_name = batch_data.get("BatchName", "")
                if batch_name != "" and not batch_name.startswith("CP"):
                    batch_names.add(batch_data["BatchName"])

    return sorted(batch_names)


def summarize_target_site_names(c2s_run_data):
    """
    Generate a list of target site names from the cells2stats report
    """
    target_site_names = set()
    for run_name in c2s_run_data:
        run_data = c2s_run_data[run_name]
        for batch_data in run_data.get("TargetStats", {}).get("Batches", []):
            for target_site_data in batch_data.get("TargetSites", []):
                target_site_name = target_site_data.get("TargetSiteName", "")
                if target_site_name != "":
                    target_site_names.add(target_site_name)
    return sorted(target_site_names)
