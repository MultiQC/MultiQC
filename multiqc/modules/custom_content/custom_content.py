"""Core MultiQC module to parse output from custom script output"""

import base64
import json
import logging
import os
import re
from collections import defaultdict
from typing import List, Dict, Union, Tuple, cast, Set, Optional, Any

import yaml

from multiqc import config, report, Plot
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, violin, heatmap, linegraph, scatter, table, box
from multiqc.plots.plotly.line import LinePlotConfig
from multiqc.plots.plotly.scatter import ScatterConfig
from multiqc.plots.plotly.box import BoxPlotConfig
from multiqc.plots.plotly.bar import BarPlotConfig
from multiqc.plots.plotly.heatmap import HeatmapConfig
from multiqc.plots.table_object import TableConfig
from multiqc.validation import ConfigValidationError

# Initialise the logger
log = logging.getLogger(__name__)


def custom_module_classes() -> List[BaseMultiqcModule]:
    """
    MultiQC Custom Content class. This module does a lot of different
    things depending on the input and is as flexible as possible.

    NB: THIS IS TOTALLY DIFFERENT TO ALL OTHER MODULES
    """

    # Dict to hold parsed data. Each key should contain a custom data type
    # e.g. output from a particular script. Note that this script may pick
    # up many different types of data from many different sources.
    # Second level keys should be 'config' and 'data'. Data key should then
    # contain sample names, and finally data.
    cust_mod_by_id: Dict[str, Dict[str, Union[Dict, List, str]]] = defaultdict(dict)

    # Dictionary to hold search patterns - start with those defined in the config
    search_patterns = ["custom_content"]

    # First - find files using patterns described in the config
    config_custom_data = getattr(config, "custom_data", {})
    mod_cust_config = {}
    for config_custom_data_id, config_custom_data_item in config_custom_data.items():
        # Check that we have a dictionary
        if not isinstance(config_custom_data_item, dict):
            log.debug(f"config.custom_data row was not a dictionary: {config_custom_data_id}")
            continue
        c_id = config_custom_data_item.get("id", config_custom_data_id)

        # Data supplied in with config (e.g. from a multiqc_config.yaml file in working directory)
        if "data" in config_custom_data_item:
            if isinstance(cust_mod_by_id.get(c_id), dict) and isinstance(cust_mod_by_id[c_id].get("data"), dict):
                d = cust_mod_by_id[c_id]["data"]
                assert isinstance(d, dict)
                d.update(config_custom_data_item["data"])
            else:
                # HTML plot type doesn't have a data sample-id key, so just take the whole chunk of data
                cust_mod_by_id[c_id]["data"] = config_custom_data_item["data"]

            cust_mod_conf = cust_mod_by_id[c_id]["config"]
            assert isinstance(cust_mod_conf, dict)
            cust_mod_conf.update({k: v for k, v in config_custom_data_item.items() if k != "data"})
            cust_mod_conf["id"] = cust_mod_conf.get("id", c_id)
            continue

        # Custom Content ID has search patterns in the config
        if report.files.get(c_id):
            config_custom_data_item["id"] = config_custom_data_item.get("id", c_id)
            cust_mod_by_id[c_id]["config"] = config_custom_data_item
            search_patterns.append(c_id)
            continue

        # Must just be configuration for a separate custom-content class
        mod_cust_config[c_id] = config_custom_data_item

    # Now go through each of the file search patterns
    bm: BaseMultiqcModule = BaseMultiqcModule(name="Custom content", anchor="custom_content")
    for config_custom_data_id in search_patterns:
        num_sp_found_files = 0
        for f in bm.find_log_files(config_custom_data_id):
            num_sp_found_files += 1

            f_extension = os.path.splitext(f["fn"])[1]

            # YAML and JSON files are the easiest
            parsed_data = None
            if f_extension == ".yaml" or f_extension == ".yml":
                try:
                    parsed_data = yaml.safe_load(f["f"])
                except Exception as e:
                    log.warning(f"Error parsing YAML file '{f['fn']}' (probably invalid YAML)")
                    log.debug(f"YAML error: {e}", exc_info=True)
                    break
                parsed_data["id"] = parsed_data.get("id", f["s_name"])
            elif f_extension == ".json":
                try:
                    parsed_data = json.loads(f["f"])
                except Exception as e:
                    log.warning(f"Error parsing JSON file '{f['fn']}' (probably invalid JSON)")
                    log.warning(f"JSON error: {e}")
                    break
                parsed_data["id"] = parsed_data.get("id", f["s_name"])
            elif f_extension == ".png" or f_extension == ".jpeg" or f_extension == ".jpg":
                image_string = base64.b64encode(f["f"].read()).decode("utf-8")
                image_format = "png" if f_extension == ".png" else "jpg"
                img_html = '<div class="mqc-custom-content-image"><img src="data:image/{};base64,{}" /></div>'.format(
                    image_format, image_string
                )
                parsed_data = {
                    "id": f["s_name"],
                    "plot_type": "image",
                    "section_name": f["s_name"].replace("_", " ").replace("-", " ").replace(".", " "),
                    "data": img_html,
                }
                cust_mod_conf = cust_mod_by_id.get(config_custom_data_id, {}).get("config", {})
                assert isinstance(cust_mod_conf, dict)
                # If the search pattern 'k' has an associated custom content section config, use it:
                parsed_data.update(cust_mod_conf)
            elif f_extension == ".html":
                parsed_data = {"id": f["s_name"], "plot_type": "html", "data": f["f"]}
                parsed_data.update(_find_html_file_header(f))

            if parsed_data is not None:
                if isinstance(parsed_data.get("data"), dict):
                    # Run sample-name cleaning on the data keys
                    parsed_data["data"] = {bm.clean_s_name(k, f): v for k, v in parsed_data["data"].items()}

                c_id = parsed_data.get("id", config_custom_data_id)
                parsed_item = parsed_data.get("data", {})
                if parsed_item:
                    cust_mod_data = cust_mod_by_id[c_id].get("data")
                    if isinstance(parsed_item, dict) and isinstance(cust_mod_data, dict):
                        cust_mod_data.update(parsed_item)
                    else:
                        cust_mod_by_id[c_id]["data"] = parsed_item
                    cust_mod_conf = cust_mod_by_id[c_id].get("config", {})
                    assert isinstance(cust_mod_conf, dict)
                    cust_mod_conf.update({j: k for j, k in parsed_data.items() if j != "data"})
                    cust_mod_by_id[c_id]["config"] = cust_mod_conf
                else:
                    log.warning(f"No data found in {f['fn']}")

            # txt, csv, tsv etc
            else:
                # Look for configuration details in the header
                m_config, non_header_lines = _find_file_header(f)
                s_name = None
                if m_config is not None:
                    c_id = m_config.get("id", config_custom_data_id)
                    # Update the base config with anything parsed from the file
                    b_config = cust_mod_by_id.get(c_id, {}).get("config", {})
                    assert isinstance(b_config, dict)
                    b_config.update(m_config)
                    # Now set the module config to the merged dict
                    m_config = dict(b_config)
                    s_name = m_config.get("sample_name")
                else:
                    c_id = config_custom_data_id
                    _m_config = cust_mod_by_id.get(c_id, {}).get("config", {})
                    assert isinstance(_m_config, dict)
                    m_config = _m_config

                # Guess sample name if not given
                if s_name is None:
                    s_name = f["s_name"]

                # Guess c_id if no information known
                if config_custom_data_id == "custom_content":
                    c_id = s_name
                    if not m_config.get("id"):
                        m_config["id"] = c_id

                # Merge with config from a MultiQC config file if we have it
                m_config.update(mod_cust_config.get(c_id, {}))

                # Add information about the file to the config dict
                if "files" not in m_config:
                    m_config["files"] = dict()
                m_config["files"].update({s_name: {"fn": f["fn"], "root": f["root"]}})

                # Guess file format if not given
                if m_config.get("file_format") is None:
                    m_config["file_format"] = _guess_file_format(f)

                # Parse data
                parsed_data, parsed_conf = _parse_txt(f, m_config, non_header_lines)
                if parsed_data is None or len(parsed_data) == 0:
                    log.warning(f"Not able to parse custom data in {f['fn']}")
                else:
                    # Did we get a new section id from the file?
                    if parsed_conf.get("id") is not None:
                        c_id = parsed_conf.get("id")
                    # heatmap - special data type
                    if isinstance(parsed_data, list):
                        cust_mod_by_id[c_id]["data"] = parsed_data
                    elif parsed_conf.get("plot_type") == "html":
                        cust_mod_by_id[c_id]["data"] = parsed_data
                    else:
                        assert isinstance(parsed_data, dict)
                        cust_data = cust_mod_by_id[c_id].get("data", {})
                        assert isinstance(cust_data, dict)
                        cust_data.update(parsed_data)
                        cust_mod_by_id[c_id]["data"] = cust_data
                    cust_config = cust_mod_by_id[c_id].get("config", {})
                    assert isinstance(cust_config, dict)
                    cust_config.update(parsed_conf)
                    cust_mod_by_id[c_id]["config"] = cust_config

        # Give log message if no files found for search pattern
        if num_sp_found_files == 0 and config_custom_data_id != "custom_content":
            log.debug(f"No samples found: custom content ({config_custom_data_id})")

    # Filter to strip out ignored sample names
    for config_custom_data_id in cust_mod_by_id:
        cust_mod_by_id[config_custom_data_id]["data"] = bm.ignore_samples(
            cust_mod_by_id[config_custom_data_id].get("data", {})
        )

    # Remove any configs that have no data
    remove_cids = [k for k in cust_mod_by_id if len(cust_mod_by_id[k].get("data", {})) == 0]
    for config_custom_data_id in remove_cids:
        del cust_mod_by_id[config_custom_data_id]

    if len(cust_mod_by_id) == 0:
        raise ModuleNoSamplesFound

    # Go through each data type
    parsed_modules: Dict[str, MultiqcModule] = dict()
    for c_id, mod_dict in cust_mod_by_id.items():
        # General Stats
        assert isinstance(mod_dict["config"], dict)
        if mod_dict["config"].get("plot_type") == "generalstats":
            assert isinstance(mod_dict["data"], dict), mod_dict["data"]
            gsheaders = mod_dict["config"].get("pconfig")
            if gsheaders is None:
                headers_set: Set[str] = set()
                for _hd in mod_dict["data"].values():
                    headers_set.update(_hd.keys())
                headers = list(headers_set)
                headers.sort()
                gsheaders = dict()
                for h in headers:
                    gsheaders[h] = dict()

            # Headers is a list of dicts
            if isinstance(gsheaders, list):
                gsheaders_dict = dict()
                for gsheader in gsheaders:
                    for col_id, col_data in gsheader.items():
                        gsheaders_dict[col_id] = col_data
                gsheaders = gsheaders_dict

            # Add namespace and description if not specified
            for m_id in gsheaders:
                if "namespace" not in gsheaders[m_id]:
                    gsheaders[m_id]["namespace"] = mod_dict["config"].get("namespace", c_id)
            log.info(f"{c_id}: Found {len(mod_dict['data'])} General Statistics columns")
            bm.general_stats_addcols(mod_dict["data"], gsheaders)

        # Initialise this new module class and append to list
        else:
            # Is this file asking to be a sub-section under a parent section?
            mod_id = mod_dict["config"].get("parent_id", mod_dict["config"].get("anchor"))
            if mod_id is None:
                mod_id = c_id
            # If we have any custom configuration from a MultiQC config file, update here
            # This is done earlier for tsv files too, but we do it here so that it overwrites what was in the file
            if mod_id in mod_cust_config:
                mod_dict["config"].update(mod_cust_config[mod_id])
            # We've not seen this module section before (normal for most custom content)
            if mod_id not in parsed_modules:
                parsed_modules[mod_id] = MultiqcModule(mod_id, mod_dict)
            else:
                # New sub-section
                parsed_modules[mod_id].update_init(mod_id, mod_dict)
            parsed_modules[mod_id].add_cc_section(c_id, mod_dict)
            if mod_dict["config"].get("plot_type") == "html":
                log.info(f"{c_id}: Found 1 sample (html)")
            elif mod_dict["config"].get("plot_type") == "image":
                log.info(f"{c_id}: Found 1 sample (image)")
            else:
                log.info(f"{c_id}: Found {len(mod_dict['data'])} samples ({mod_dict['config'].get('plot_type')})")

    # Sort sections if we have a config option for order
    mod_order = getattr(config, "custom_content", {}).get("order", [])
    # after each element, also add a version with a "-module" next to it for back-compat with < 1.24
    mod_order = [x for x in mod_order for x in [x, re.sub("-module$", "", x), f"{x}-module"]]
    modules__not_in_order: List[BaseMultiqcModule] = [
        parsed_mod for parsed_mod in parsed_modules.values() if parsed_mod.anchor not in mod_order
    ]
    modules_in_order = [
        parsed_mod for mod_id in mod_order for parsed_mod in parsed_modules.values() if parsed_mod.anchor == mod_id
    ]
    sorted_modules = modules_in_order + modules__not_in_order

    # If we only have General Stats columns then there are no module outputs
    if len(sorted_modules) == 0:
        cfgs: List[Dict] = []
        for cust_mod_v in cust_mod_by_id.values():
            assert isinstance(cust_mod_v["config"], dict)
            cfgs.append(cust_mod_v["config"])
        if len(cust_mod_by_id) >= 1 and all(cfg.get("plot_type") == "generalstats" for cfg in cfgs):
            sorted_modules = [bm]
        else:
            raise ModuleNoSamplesFound

    return sorted_modules


class MultiqcModule(BaseMultiqcModule):
    """Module class, used for each custom content type"""

    def __init__(self, c_id: str, mod: Dict[str, Union[Dict, List, str]]):
        modname = c_id.replace("_", " ").title()
        assert isinstance(mod["config"], dict)
        mod_info = mod["config"].get("description")
        if "parent_name" in mod["config"]:
            assert isinstance(mod["config"]["parent_name"], str)
            modname = mod["config"]["parent_name"]
            mod_info = mod["config"].get("parent_description")
        elif "section_name" in mod["config"]:
            assert isinstance(mod["config"]["section_name"], str)
            modname = mod["config"]["section_name"]
        if modname == "" or modname is None:
            modname = "Custom Content"

        anchor = mod["config"].get("section_anchor", c_id)

        super(MultiqcModule, self).__init__(
            name=modname,
            anchor=anchor,
            href=mod["config"].get("section_href"),
            info=mod_info,
            extra=mod["config"].get("extra"),
            doi=mod["config"].get("doi"),
        )

        if "custom_content" in config.run_modules:
            # To allow file_search.include_or_exclude_modules() correctly filter these modules
            config.custom_content_modules.append(anchor)

    def update_init(self, c_id, mod):
        """
        This function runs when we have already initialised a module
        and this is a subsequent file that will be another subsection.

        So most info should already be initialised properly.
        We check if anything is empty and if we have it set in this file
        we set it here.
        """

        if self.info is None or self.info == "":
            self.info = mod["config"].get("parent_description")
        if self.extra is None or self.info == "":
            self.extra = mod["config"].get("extra", None)
        # This needs overwriting again as it has already run on init
        self.intro = self._get_intro()

    def add_cc_section(self, c_id, mod):
        section_name = mod["config"].get("section_name", c_id.replace("_", " ").title())
        if section_name == "" or section_name is None:
            section_name = "Custom Content"

        if self.id == c_id:  # make sure the anchors are unique
            c_id = c_id + "-section"

        pconfig = mod["config"].get("pconfig", {})
        if pconfig.get("id") is None:
            pconfig["id"] = f"{c_id}-plot"
        if pconfig["id"] == c_id or pconfig["id"] == mod["config"].get("id"):
            pconfig["id"] += "-plot"
        if pconfig.get("title") is None:
            pconfig["title"] = section_name

        plot: Optional[Union[Plot, str]] = None
        content = None

        plot_type = mod["config"].get("plot_type")

        # Heatmap
        if plot_type == "heatmap":
            plot = heatmap.plot(
                mod["data"], mod["config"].get("xcats"), mod["config"].get("ycats"), pconfig=HeatmapConfig(**pconfig)
            )
        else:
            if not isinstance(mod["data"], list):
                mod["data"] = [mod["data"]]

            # Try to coerce x-axis to numeric
            if plot_type in ["linegraph", "scatter"]:
                try:
                    mod["data"] = [{k: {float(x): v[x] for x in v} for k, v in ds.items()} for ds in mod["data"]]
                except ValueError:
                    pass

            # Table
            if plot_type == "table":
                headers = mod["config"].get("headers")

                # handle some legacy fields for backwards compat
                sort_rows = pconfig.pop("sortRows", None)
                if sort_rows is not None:
                    pconfig["sort_rows"] = sort_rows
                no_violin = pconfig.pop("no_beeswarm", None)
                if no_violin is not None:
                    pconfig["no_violin"] = no_violin

                plot = table.plot(mod["data"], headers=headers, pconfig=pconfig)

            # Bar plot
            elif plot_type == "bargraph":
                mod["data"] = [{str(k): v for k, v in ds.items()} for ds in mod["data"]]
                plot = bargraph.plot(mod["data"], mod["config"].get("categories"), pconfig=BarPlotConfig(**pconfig))

            # Line plot
            elif plot_type == "linegraph":
                plot = linegraph.plot(mod["data"], pconfig=LinePlotConfig(**pconfig))

            # Scatter plot
            elif plot_type == "scatter":
                plot = scatter.plot(mod["data"], pconfig=ScatterConfig(**pconfig))

            # Box plot
            elif plot_type == "box":
                plot = box.plot(mod["data"], pconfig=BoxPlotConfig(**pconfig))

            # Violin plot
            elif plot_type in ["violin", "beeswarm"]:
                plot = violin.plot(mod["data"], pconfig=TableConfig(**pconfig))

            # Raw HTML
            elif plot_type == "html":
                if len(mod["data"]) > 1:
                    log.warning(f"HTML plot type found with more than one dataset in {c_id}")
                content = mod["data"][0]

            # Raw image file as html
            elif plot_type == "image":
                if len(mod["data"]) > 1:
                    log.warning(f"Image plot type found with more than one dataset in {c_id}")
                content = mod["data"][0]

            # Not supplied
            elif plot_type is None:
                log.warning(f"Plot type not found for content ID '{c_id}'")

            # Not recognised
            else:
                log.warning(
                    "Error - custom content plot type '{}' not recognised for content ID {}".format(
                        mod["config"].get("plot_type"), c_id
                    )
                )

        if plot is not None:
            for i, ds in enumerate(mod["data"]):
                # Save the data if it's not a html string
                if not isinstance(plot, str):
                    did = plot.pconfig.id
                    if i > 0:
                        did = f"{did}_{i}"
                    self.write_data_file(ds, f"multiqc_{did}")
                    plot.pconfig.save_data_file = False

        # Don't use exactly the same title / description text as the main module
        if section_name == self.name:
            section_name = None

        self.add_section(name=section_name, anchor=c_id, plot=plot, content=content)


def _find_file_header(f) -> Tuple[Optional[Dict], List[str]]:
    # Collect commented out header lines
    hlines = []
    other_lines = []
    for line in f["f"].splitlines():
        if line.startswith("#"):
            hlines.append(line[1:])
        else:
            other_lines.append(line)

    # Check if the last header line is the '#'-commented column names
    sep = None
    if f["fn"].endswith(".tsv"):
        sep = "\t"
    elif f["fn"].endswith(".csv"):
        sep = ","
    if (
        sep
        and len(hlines) > 0
        and len(other_lines) > 0
        and ":" not in hlines[-1]
        and len(hlines[-1].split(sep)) == len(other_lines[0].split(sep))
    ):
        other_lines = [hlines.pop()] + other_lines

    if len(hlines) == 0:
        return None, other_lines

    try:
        hconfig = yaml.safe_load("\n".join(hlines))
    except yaml.YAMLError:
        raise ConfigValidationError(
            f"Could not parse comment file header for MultiQC custom content: {f['fn']}. "
            + "Note that everything behind a comment character '#' is expected to in YAML format. "
            + "To provide column names in a TSV or CSV file, put them as the first raw without fencing it with a '#'.",
            "custom_content",
        )
    else:
        if not isinstance(hconfig, dict):
            raise ConfigValidationError(
                "Custom Content comment file header looked wrong. It's expected to "
                + f"be parsed to a dict, got {type(hconfig)}: {hconfig}",
                "custom_content",
            )
    return hconfig, other_lines


def _find_html_file_header(f):
    """Look for a HTML comment config at the start of a custom content HTML file"""
    if f["f"].lstrip().startswith("<!--"):
        match = re.search(r"^\<\!\-\-((?:.|\n|\r)*?)-->", f["f"].lstrip())
        if match:
            comment = match.group(1)
            if comment:
                try:
                    return yaml.load(comment, Loader=yaml.SafeLoader)
                except Exception as e:
                    log.debug(f"Found Custom Content HTML comment, but couldn't load as YAML: {e}", exc_info=True)
                    log.debug(f"Comment:\n{comment}")
    return {}


def _guess_file_format(f):
    """
    Tries to guess file format, first based on file extension (csv / tsv),
    then by looking for common column separators in the first 10 non-commented lines.
    Splits by tab / comma / space and counts resulting number of columns. Finds the most
    common column count, then compared how many lines had this number.
    e.g. if tab, all 10 lines should have x columns when split by tab.
    Returns: csv | tsv | spaces   (spaces by default if all else fails)
    """
    filename, file_extension = os.path.splitext(f["fn"])
    tabs = []
    commas = []
    spaces = []
    j = 0
    for line in f["f"].splitlines():
        if not line.startswith("#"):
            j += 1
            tabs.append(len(line.split("\t")))
            commas.append(len(line.split(",")))
            spaces.append(len(line.split()))
        if j == 10:
            break
    tab_mode = max(set(tabs), key=tabs.count)
    commas_mode = max(set(commas), key=commas.count)
    spaces_mode = max(set(spaces), key=spaces.count)
    tab_lc = tabs.count(tab_mode) if tab_mode > 1 else 0
    commas_lc = commas.count(commas_mode) if commas_mode > 1 else 0
    spaces_lc = spaces.count(spaces_mode) if spaces_mode > 1 else 0
    if tab_lc == j:
        return "tsv"
    elif commas_lc == j:
        return "csv"
    else:
        if tab_lc > commas_lc and tab_lc > spaces_lc:
            return "tsv"
        elif commas_lc > tab_lc and commas_lc > spaces_lc:
            return "csv"
        elif spaces_lc > tab_lc and spaces_lc > commas_lc:
            return "spaces"
        else:
            if tab_mode == commas_lc and tab_mode > spaces_lc:
                if tab_mode > commas_mode:
                    return "tsv"
                else:
                    return "csv"
    return "spaces"


def _parse_txt(f, conf: Dict, non_header_lines: List[str]) -> Tuple[Union[str, Dict, List, None], Dict]:
    # Split the data into a list of lists by column
    sep = None
    if conf["file_format"] == "csv":
        sep = ","
    if conf["file_format"] == "tsv":
        sep = "\t"

    # Check for special case - HTML
    if conf.get("plot_type") == "html":
        lines: List[str] = []
        for line in non_header_lines:
            if line:
                lines.append(line)
        return "\n".join(lines), conf

    # Not HTML, need to parse data
    matrix: List[List[Union[str, float, int]]] = []
    ncols = None
    for line in non_header_lines:
        if line:
            sections = cast(List[Any], line.split(sep))
            matrix.append(sections)
            if ncols is None:
                ncols = len(sections)
            elif ncols != len(sections):
                log.warning(f"Inconsistent number of columns found in {f['fn']}! Skipping..")
                return None, conf

    # Convert values to floats if we can
    first_row_str = 0
    for i, sections in enumerate(matrix):
        for j, v in enumerate(sections):
            try:
                v = float(v)
            except ValueError:
                pass
            if isinstance(v, str):
                if (v.startswith('"') and v.endswith('"')) or (v.startswith("'") and v.endswith("'")):
                    v = v[1:-1]
                # Count strings in first row (header?)
                if i == 0:
                    first_row_str += 1
            matrix[i][j] = v

    all_numeric = all(isinstance(v, (float, int)) for s in matrix[1:] for v in s[1:])

    # General stat info files - expected to have at least 2 rows (first row always being the header)
    # and have at least 2 columns (first column always being sample name)
    if conf.get("plot_type") == "generalstats" and len(matrix) >= 2 and ncols and ncols >= 2:
        data_ddict: Dict[str, Dict] = defaultdict(dict)
        for i, section in enumerate(matrix[1:], 1):
            for j, v in enumerate(section[1:], 1):
                s_name = str(section[0])
                data_ddict[s_name][matrix[0][j]] = v
        return data_ddict, conf

    # Heatmap: Number of headers == number of lines
    if conf.get("plot_type") is None and first_row_str == len(non_header_lines) and all_numeric:
        conf["plot_type"] = "heatmap"
    if conf.get("plot_type") == "heatmap":
        conf["xcats"] = matrix[0][1:]
        conf["ycats"] = [s[0] for s in matrix[1:]]
        data_list: List = [s[1:] for s in matrix[1:]]
        return data_list, conf

    # Header row of strings, or configured as table
    if first_row_str == len(matrix[0]) or conf.get("plot_type") == "table":
        data_ddict = dict()
        for s in matrix[1:]:
            s_name = str(s[0])
            data_ddict[s_name] = dict()
            for i, v in enumerate(s[1:]):
                cat = str(matrix[0][i + 1])
                data_ddict[s_name][cat] = v
        # Bar graph or table - if numeric data, go for bar graph
        if conf.get("plot_type") is None:
            allfloats = True
            for r in matrix[1:]:
                for v in r[1:]:
                    allfloats = allfloats and isinstance(v, float)
            if allfloats:
                conf["plot_type"] = "bargraph"
            else:
                conf["plot_type"] = "table"
        # Set table col_1 header
        col_name = str(matrix[0][0])
        if conf.get("plot_type") == "table" and col_name.strip() != "":
            conf["pconfig"] = conf.get("pconfig", {})
            if not conf["pconfig"].get("col1_header"):
                conf["pconfig"]["col1_header"] = col_name.strip()
        # Return parsed data
        if conf.get("plot_type") == "bargraph" or conf.get("plot_type") == "table":
            return data_ddict, conf
        else:
            data_ddict = dict()  # reset

    # Scatter plot: First row is  str : num : num
    if (
        conf.get("plot_type") is None
        and len(matrix[0]) == 3
        and not isinstance(matrix[0][0], float)
        and isinstance(matrix[0][1], float)
        and isinstance(matrix[0][2], float)
    ):
        conf["plot_type"] = "scatter"

    if conf.get("plot_type") == "scatter":
        dicts: Dict[str, Dict[str, float]] = dict()
        for s in matrix:
            try:
                dicts[str(s[0])] = {"x": float(s[1]), "y": float(s[2])}
            except (IndexError, ValueError):
                pass
        return dicts, conf

    # Single sample line / bar graph - first row has two columns
    if len(matrix[0]) == 2:
        # Line graph - num : num
        if conf.get("plot_type") is None and isinstance(matrix[0][0], float) and isinstance(matrix[0][1], float):
            conf["plot_type"] = "linegraph"
        # Bar graph - str : num
        if conf.get("plot_type") is None and not isinstance(matrix[0][0], float) and isinstance(matrix[0][1], float):
            conf["plot_type"] = "bargraph"

        # Data structure is the same
        if conf.get("plot_type") == "linegraph" or conf.get("plot_type") == "bargraph":
            # Set section id based on directory if not known
            if conf.get("id") is None:
                conf["id"] = os.path.basename(f["root"])
            data_dict: Dict = dict()
            for s in matrix:
                data_dict[s[0]] = s[1]
            return {f["s_name"]: data_dict}, conf

    # Multi-sample line graph: No header row, str : lots of num columns
    if conf.get("plot_type") is None and len(matrix[0]) > 4 and all_numeric:
        conf["plot_type"] = "linegraph"

    if conf.get("plot_type") == "linegraph":
        data_ddict = dict()
        # If the first row has no header, use it as axis labels
        x_labels = []
        s_name = str(matrix[0][0])
        if s_name.strip() == "":
            x_labels = matrix.pop(0)[1:]
        # Use 1..n range for x values
        for s in matrix:
            name = str(s[0])
            data_ddict[name] = dict()
            for i, v in enumerate(s[1:]):
                try:
                    x_val = x_labels[i]
                    try:
                        x_val = float(x_val)
                    except ValueError:
                        pass
                except IndexError:
                    x_val = i + 1
                data_ddict[name][x_val] = v
        return data_ddict, conf

    # Got to the end and haven't returned. It's a mystery, capn'!
    log.debug(
        f"Not able to figure out a plot type for '{f['fn']}' "
        + "plot type = {}, all numeric = {}, first row str = {}".format(
            conf.get("plot_type"), all_numeric, first_row_str
        )
    )
    return None, conf
