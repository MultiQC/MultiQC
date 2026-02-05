import logging
from typing import Dict, List, Tuple, Optional

from multiqc import config

log = logging.getLogger(__name__)


def table_data_and_headers(
    rows_list: List[Tuple[str, str]],
    help: List[Tuple[str, List[str]]],
    namespace: Optional[str] = None,
) -> Tuple[Dict, Dict]:
    """Update the data dict and headers dict

    Args:
        rows_list: List of tuples containing column names and data values
        help: List of tuples containing help text for each column
        namespace: Optional namespace prefix (e.g., "ATAC", "GEX") to prevent key collisions.
                   Keys in the returned dicts will be prefixed with "{namespace}_" if provided,
                   while display titles remain clean (without prefix).

    Returns:
        Tuple of (table dict, headers dict) with optionally prefixed keys
    """
    table = dict()
    headers = dict()

    help_dict = {lst[0]: lst[1][0].strip(".") for lst in help}
    for col_name, col_data in rows_list:
        # Sanitize numeric data
        is_percentage = "%" in col_data
        col_data = col_data.replace(",", "").replace("%", "")
        is_integer = col_data.isdigit()

        # Convert to float when possible
        try:
            col_data = float(col_data)
        except ValueError:
            col_data = col_data

        # Apply namespace prefix to key if provided
        key = f"{namespace}_{col_name}" if namespace else col_name

        table[key] = col_data

        # Assign shared/regular keys
        if col_name == "Sequenced read pairs":
            headers[key] = {
                "title": col_name,
                "description": f"{help_dict[col_name]} ({config.read_count_desc})",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "format": "{:,.0f}",
            }
        else:
            headers[key] = {
                "title": col_name,
                "description": help_dict[col_name],
            }
        if col_name == "Estimated number of cells":
            headers[key]["shared_key"] = "cell_count"
            headers[key]["title"] = "Est. cells"
            headers[key]["description"] += " (which estimates the number of cells)"

        if col_name == "Fraction of high-quality fragments in cells":
            headers[key]["title"] = "High-qual fragments"

        if is_integer:
            headers[key]["format"] = "{:,.0f}"

        if is_percentage:
            headers[key].update(
                {
                    "suffix": "%",
                    "max": 100,
                    "min": 0,
                }
            )

    return table, headers


def set_hidden_cols(headers, col_names):
    """Set the hidden columns

    Note: This function handles both prefixed keys (e.g., "ATAC_Percent duplicates")
    and unprefixed keys. It will try to find the key in headers and set it as hidden.
    """

    for col_name in col_names:
        try:
            headers[col_name]["hidden"] = True
        except KeyError:
            pass

    return headers


def subset_header(data, cols, namespace=None):
    """Subsets the headers to only columns in cols. Adds colour and namespace

    Args:
        data: The complete headers dictionary (may contain prefixed keys)
        cols: Dictionary mapping column names to color scales
        namespace: Optional namespace for the plot display

    Note: This function handles prefixed keys in the data dict. It looks up keys
    from cols in the data dict, which may have namespace prefixes applied.
    """

    headers = dict()
    for key, val in cols.items():
        # Key lookup: try exact match first, then check if it exists in data
        if key in data:
            headers[key] = data[key]
        else:
            # If exact match fails, skip this key (it may not exist in this dataset)
            log.debug(f"Key '{key}' not found in headers data")
            continue

        headers[key]["scale"] = val
        headers[key]["namespace"] = namespace

    return headers


def extract_plot_data(data):
    """Extracts plot data"""

    plot_data = dict()
    axes_data = data["data"][0]
    plot_data = dict(zip(axes_data["x"], axes_data["y"]))
    return plot_data
