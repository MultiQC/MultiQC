def clean_title_case(col_id):
    title = col_id.title() if col_id[0:1].islower() else col_id
    for _str in ["Bc", "bc", "Umi", "Igk", "Igh", "Igl", "Vj", "q30"]:
        title = title.replace(_str, _str.upper())
    return title


def update_data_and_headers(data, headers, new_data, new_headers, colors, prefix, int_cols=()):
    """update the data dict and headers dict

    :param: int_cols columns to be shown as integers in the table
    """

    for col_name, col_data in new_data:
        if col_name in new_headers:
            # Sanitize numeric data
            is_percentage = "%" in col_data
            col_data = col_data.replace(",", "").replace("%", "")

            # Convert to float when possible
            if col_data == "None":
                col_data = None
            else:
                try:
                    col_data = float(col_data)
                except ValueError:
                    pass

            col_id = new_headers[col_name]
            data[col_id] = col_data
            headers[col_id] = {
                "rid": "{}_{}".format(prefix, col_id.replace(" ", "_").replace("/", "_")),
                "title": clean_title_case(col_id),
                "description": col_name,
                "namespace": f"Space Ranger {prefix}",
                "scale": colors.get(col_id, "RdYlGn" if is_percentage else "GnBu"),
            }
            if is_percentage:
                headers[col_id].update(
                    {
                        "suffix": "%",
                        "max": 100,
                        "min": 0,
                    }
                )

            if col_id in int_cols:
                headers[col_id]["format"] = "{:,.0f}"


def set_hidden_cols(headers, col_names):
    """Set the hidden columns"""

    for col_name in col_names:
        try:
            headers[col_name]["hidden"] = True
        except KeyError:
            pass

    return headers


def parse_bcknee_data(data, s_name, max_idx=1000):
    """parse data for bc knee plot from spaceranger dict"""

    value_dict = dict()
    for idx, data_series in enumerate(data):
        if idx > max_idx:
            break
        if len(data_series["x"]) == 0:
            continue
        id = f"{s_name}_{data_series['name']}"
        if id not in value_dict.keys():
            value_dict[id] = dict()
        value_dict[id].update(transform_data(data_series))

    return value_dict


def transform_data(data):
    """Transform x:list,y:list data to a dict of x_val:y_val"""

    return {x: y for x, y in zip(data["x"], data["y"])}
