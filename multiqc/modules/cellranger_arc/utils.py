from multiqc import config


def update_dict(rows_list, help):
    """update the data dict and headers dict"""
    table = dict()
    headers = dict()

    help_dict = {lst[0]: lst[1][0] for lst in help}
    for col_name, col_data in rows_list:
        # Sanitize numeric data
        is_percentage = "%" in col_data
        col_data = col_data.replace(",", "").replace("%", "")

        # Convert to float when possible
        try:
            col_data = float(col_data)
        except ValueError:
            col_data = col_data

        table[col_name] = col_data

        # Assign shared/regular keys
        if col_name == "Sequenced read pairs":
            headers[col_name] = {
                "title": f"{config.read_count_prefix} {col_name}",
                "description": f"{help_dict[col_name]} ({config.read_count_desc})",
                "modify": lambda x: x * config.read_count_multiplier,
                "shared_key": "read_count",
                "format": "{:,.0f}",
            }
        else:
            headers[col_name] = {
                "title": col_name,
                "description": help_dict[col_name],
            }
        if col_name == "Estimated number of cells":
            headers[col_name]["shared_key"] = "cell_count"
            headers[col_name]["format"] = "{:,.0f}"

        if is_percentage:
            headers[col_name].update(
                {
                    "suffix": "%",
                    "max": 100,
                    "min": 0,
                }
            )

    return table, headers


def set_hidden_cols(headers, col_names):
    """Set the hidden columns"""

    for col_name in col_names:
        try:
            headers[col_name]["hidden"] = True
        except KeyError:
            pass

    return headers


def subset_header(data, cols, namespace=None):
    """Subsets the headers to only columns in cols. Adds colour and namespace"""

    headers = dict()
    for key, val in cols.items():
        headers[key] = data[key]
        headers[key]["scale"] = val
        headers[key]["namespace"] = namespace

    return headers


def extract_plot_data(data):
    """Extracts plot data"""

    plot_data = dict()
    axes_data = data["data"][0]
    plot_data = dict(zip(axes_data["x"], axes_data["y"]))
    return plot_data