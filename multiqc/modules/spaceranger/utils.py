from typing import Dict, Union


def clean_title_case(col_id):
    title = col_id.title() if col_id[0:1].islower() else col_id
    for _str in ["Bc", "bc", "Umi", "Igk", "Igh", "Igl", "Vj", "q30"]:
        title = title.replace(_str, _str.upper())
    return title


def populate_data_and_headers(
    headers_to_update: Dict[str, Dict[str, Union[int, str, float, None]]],
    new_data,
    new_headers,
    colors,
    prefix,
    int_cols=(),
) -> Dict[str, Union[int, str, float, None]]:
    """
    Populate the data dict and headers dict.

    `int_cols` columns to be shown as integers in the table
    """

    val_by_metric = dict()

    for col_name, value in new_data:
        if col_name in new_headers:
            # Sanitize numeric data
            is_percentage = "%" in value
            value = value.replace(",", "").replace("%", "")

            # Convert to float when possible
            if value == "None":
                value = None
            else:
                try:
                    value = float(value)
                except ValueError:
                    pass

            col_id = new_headers[col_name]
            val_by_metric[col_id] = value
            if col_id not in headers_to_update:
                headers_to_update[col_id] = {
                    "rid": "{}_{}".format(prefix, col_id.replace(" ", "_").replace("/", "_")),
                    "title": clean_title_case(col_id),
                    "description": col_name,
                    "namespace": f"Space Ranger {prefix}",
                    "scale": colors.get(col_id, "RdYlGn" if is_percentage else "GnBu"),
                }
                if is_percentage:
                    headers_to_update[col_id].update(
                        {
                            "suffix": "%",
                            "max": 100,
                            "min": 0,
                        }
                    )
                if col_id in int_cols:
                    headers_to_update[col_id]["format"] = "{:,.0f}"

    return val_by_metric


def set_hidden_cols(headers, col_names):
    """Set the hidden columns"""

    for col_name in col_names:
        try:
            headers[col_name]["hidden"] = True
        except KeyError:
            pass

    return headers


def transform_data(data: Dict[str, list]) -> Dict[str, Union[int, str, float]]:
    """Transform x:list,y:list data to a dict of x_val:y_val"""

    return {x: y for x, y in zip(data["x"], data["y"])}
