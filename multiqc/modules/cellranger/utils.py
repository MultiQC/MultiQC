import logging
from typing import Dict, Union, List, Tuple, Any


class SectionData:
    data: Dict[str, Dict]
    headers: Dict[str, Any]

    def __init__(self):
        self.data = dict()
        self.headers = dict()

    def update(self, data: Dict[str, Dict], headers: Dict[str, Dict]):
        self.data.update(data)
        self.headers.update(headers)

    def update_sample(self, data: Dict[str, Any], headers: Dict[str, Dict], sample: str):
        if sample not in self.data:
            self.data[sample] = {}
        self.data[sample].update(data)
        self.headers.update(headers)


def clean_title_case(col_id):
    title = col_id.title() if col_id[0:1].islower() else col_id
    for str in ["Bc", "bc", "Umi", "Igk", "Igh", "Igl", "Vj", "q30"]:
        title = title.replace(str, str.upper())
    return title


def update_dict(table, headers, rows_list, col_map, colours, namespace):
    """update the data dict and headers dict"""

    for col_name, col_data in rows_list:
        if col_name in col_map:
            # Sanitize numeric data
            is_percentage = "%" in col_data
            col_data = col_data.replace(",", "").replace("%", "")

            # Convert to float when possible
            try:
                col_data = float(col_data)
            except ValueError:
                col_data = col_data

            col_id = col_map[col_name]
            table[col_id] = col_data
            headers[col_id] = {
                "rid": f"{namespace}_{col_id.replace(' ', '_').replace('/', '_')}".lower(),
                "title": clean_title_case(col_id),
                "description": col_name,
                "namespace": namespace,
                "scale": colours.get(col_id, "RdYlGn" if is_percentage else "GnBu"),
            }
            if is_percentage:
                headers[col_id].update(
                    {
                        "suffix": "%",
                        "max": 100,
                        "min": 0,
                    }
                )

            # Assign shared keys
            if col_id == "estimated cells":
                headers[col_id]["shared_key"] = "cell_count"
                headers[col_id]["format"] = "{:,.0f}"
            if col_id.endswith("reads/cell"):
                headers[col_id]["shared_key"] = "reads_per_cell"
                headers[col_id]["format"] = "{:,.0f}"

    return table, headers


def set_hidden_cols(headers, col_names):
    """Set the hidden columns"""

    for col_name in col_names:
        try:
            headers[col_name]["hidden"] = True
        except KeyError:
            pass

    return headers


def parse_bcknee_data(data, s_name, max_idx=1000) -> Dict[str, Dict[str, Union[int, str, float]]]:
    """parse data for bc knee plot from dict"""

    value_dict: Dict[str, Dict[str, Union[int, str, float]]] = dict()
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


def transform_data(data: Dict[str, List]) -> Dict[str, Union[int, str, float]]:
    """Transform x:list,y:list data to a dict of x_val:y_val"""

    value_dict = dict()
    for idx, row in enumerate(data["x"]):
        if row > 0 and data["y"][idx] > 0:
            value_dict[row] = data["y"][idx]

    return value_dict


def update_sections_tuple(original_tuple: Tuple, delta_tuple: Tuple):
    """Helper to apply Section.update(delta.data, delta.headers) for each tuple element"""
    if len(delta_tuple) != len(original_tuple):
        logging.critical("len(delta_tuple) != len(original_tuple), Likely to be a bug")

    for original, delta in zip(original_tuple, delta_tuple):
        original.update(delta.data, delta.headers)


def resolve_dict(data: dict, path: list):
    """Return the contents of a dictionary after traversing the path provided"""

    cdata = data
    for entry in path:
        cdata = cdata[entry]

    return cdata


def update_data_and_header(data: Dict, header: Dict, color_dict: Dict):
    """Transform data and headers of a table"""

    for key in data:
        new_header = {"scale": color_dict.get(key, "GnBu")}
        value = data[key].replace(",", "")
        if "%" in value:
            new_header.update({"suffix": "%", "max": 100, "min": 0})
            value = value.replace("%", "").strip()
        try:
            value = float(value)
            value = int(value) if value % 1 == 0 else value
        except ValueError:
            pass
        data[key] = value
        header[key].update(new_header)
    return data, header


def combine_data_and_headers(dictionaries: List[Tuple]):
    """Method to combine multiple data dictionaries and headers into one"""

    ret_data: Dict[str, Dict] = {}
    ret_headers: Dict[str, Any] = {}
    for data, headers, cols in dictionaries:
        for sample in data:
            if sample not in ret_data:
                ret_data[sample] = {}
            for col in cols:
                if col in data[sample]:
                    ret_data[sample][col] = data[sample][col]
                    ret_headers[col] = headers[col]

    return ret_data, ret_headers
