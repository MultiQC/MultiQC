""" Core MultiQC module to parse output from custom script output """


import base64
import json
import logging
import os
import re
from collections import defaultdict
from typing import List, Dict

import yaml

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, violin, heatmap, linegraph, scatter, table
from multiqc.utils import report

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
    cust_mod_by_id: Dict[str, Dict[str, Dict]] = defaultdict(lambda: defaultdict(lambda: dict()))

    # Dictionary to hold search patterns - start with those defined in the config
    search_patterns = ["custom_content"]

    # First - find files using patterns described in the config
    config_data = getattr(config, "custom_data", {})
    mod_cust_config = {}
    for k, f in config_data.items():
        # Check that we have a dictionary
        if not isinstance(f, dict):
            log.debug(f"config.custom_data row was not a dictionary: {k}")
            continue
        c_id = f.get("id", k)

        # Data supplied in with config (e.g. from a multiqc_config.yaml file in working directory)
        if "data" in f:
            try:
                cust_mod_by_id[c_id]["data"].update(f["data"])
            except ValueError:
                # HTML plot type doesn't have a data sample-id key, so just take the whole chunk of data
                cust_mod_by_id[c_id]["data"] = f["data"]
            cust_mod_by_id[c_id]["config"].update({k: v for k, v in f.items() if k != "data"})
            cust_mod_by_id[c_id]["config"]["id"] = cust_mod_by_id[c_id]["config"].get("id", c_id)
            continue

        # Custom Content ID has search patterns in the config
        if c_id in report.files:
            cust_mod_by_id[c_id]["config"] = f
            cust_mod_by_id[c_id]["config"]["id"] = cust_mod_by_id[c_id]["config"].get("id", c_id)
            search_patterns.append(c_id)
            continue

        # Must just be configuration for a separate custom-content class
        mod_cust_config[c_id] = f

    # Now go through each of the file search patterns
    bm = BaseMultiqcModule(name="Custom content", anchor="custom_content")
    for k in search_patterns:
        num_sp_found_files = 0
        for f in bm.find_log_files(k):
            num_sp_found_files += 1
            # Handle any exception without messing up for remaining custom content files
            try:
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
                    img_html = (
                        '<div class="mqc-custom-content-image"><img src="data:image/{};base64,{}" /></div>'.format(
                            image_format, image_string
                        )
                    )
                    parsed_data = {
                        "id": f["s_name"],
                        "plot_type": "image",
                        "section_name": f["s_name"].replace("_", " ").replace("-", " ").replace(".", " "),
                        "data": img_html,
                    }
                    # If the search pattern 'k' has an associated custom content section config, use it
                    parsed_data.update(cust_mod_by_id.get(k, {}).get("config", {}))
                elif f_extension == ".html":
                    parsed_data = {"id": f["s_name"], "plot_type": "html", "data": f["f"]}
                    parsed_data.update(_find_html_file_header(f))

                if parsed_data is not None:
                    if isinstance(parsed_data.get("data"), dict):
                        # Run sample-name cleaning on the data keys
                        parsed_data["data"] = {bm.clean_s_name(k, f): v for k, v in parsed_data["data"].items()}

                    c_id = parsed_data.get("id", k)
                    if len(parsed_data.get("data", {})) > 0:
                        if isinstance(parsed_data["data"], dict):
                            cust_mod_by_id[c_id]["data"].update(parsed_data["data"])
                        else:
                            cust_mod_by_id[c_id]["data"] = parsed_data["data"]
                        cust_mod_by_id[c_id]["config"].update({j: k for j, k in parsed_data.items() if j != "data"})
                    else:
                        log.warning(f"No data found in {f['fn']}")

                # txt, csv, tsv etc
                else:
                    # Look for configuration details in the header
                    m_config = _find_file_header(f)
                    s_name = None
                    if m_config is not None:
                        c_id = m_config.get("id", k)
                        # Update the base config with anything parsed from the file
                        b_config = cust_mod_by_id.get(c_id, {}).get("config", {})
                        b_config.update(m_config)
                        # Now set the module config to the merged dict
                        m_config = dict(b_config)
                        s_name = m_config.get("sample_name")
                    else:
                        c_id = k
                        m_config = cust_mod_by_id.get(c_id, {}).get("config", {})

                    # Guess sample name if not given
                    if s_name is None:
                        s_name = f["s_name"]

                    # Guess c_id if no information known
                    if k == "custom_content":
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
                    try:
                        parsed_data, conf = _parse_txt(f, m_config)
                        if parsed_data is None or len(parsed_data) == 0:
                            log.warning(f"Not able to parse custom data in {f['fn']}")
                        else:
                            # Did we get a new section id from the file?
                            if conf.get("id") is not None:
                                c_id = conf.get("id")
                            # heatmap - special data type
                            if isinstance(parsed_data, list):
                                cust_mod_by_id[c_id]["data"] = parsed_data
                            elif conf.get("plot_type") == "html":
                                cust_mod_by_id[c_id]["data"] = parsed_data
                            else:
                                cust_mod_by_id[c_id]["data"].update(parsed_data)
                            cust_mod_by_id[c_id]["config"].update(conf)
                    except (IndexError, AttributeError, TypeError):
                        log.error(f"Unexpected parsing error for {f['fn']}", exc_info=True)
                        raise  # testing
            except Exception as e:
                log.error(f"Uncaught exception raised for file '{f['fn']}'")
                log.exception(e)

        # Give log message if no files found for search pattern
        if num_sp_found_files == 0 and k != "custom_content":
            log.debug(f"No samples found: custom content ({k})")

    # Filter to strip out ignored sample names
    for k in cust_mod_by_id:
        cust_mod_by_id[k]["data"] = bm.ignore_samples(cust_mod_by_id[k]["data"])

    # Remove any configs that have no data
    remove_cids = [k for k in cust_mod_by_id if len(cust_mod_by_id[k]["data"]) == 0]
    for k in remove_cids:
        del cust_mod_by_id[k]

    if len(cust_mod_by_id) == 0:
        raise ModuleNoSamplesFound

    # Go through each data type
    parsed_modules = dict()
    for c_id, mod in cust_mod_by_id.items():
        # General Stats
        if mod["config"].get("plot_type") == "generalstats":
            gsheaders = mod["config"].get("pconfig")
            if gsheaders is None:
                headers = set()
                for d in mod["data"].values():
                    headers.update(d.keys())
                headers = list(headers)
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
                    gsheaders[m_id]["namespace"] = mod["config"].get("namespace", c_id)
            log.info(f"{c_id}: Found {len(mod['data'])} General Statistics columns")
            bm.general_stats_addcols(mod["data"], gsheaders)

        # Initialise this new module class and append to list
        else:
            # Is this file asking to be a sub-section under a parent section?
            mod_id = mod["config"].get("parent_id", c_id)
            # If we have any custom configuration from a MultiQC config file, update here
            # This is done earlier for tsv files too, but we do it here so that it overwrites what was in the file
            if mod_id in mod_cust_config:
                mod["config"].update(mod_cust_config[mod_id])
            # We've not seen this module section before (normal for most custom content)
            if mod_id not in parsed_modules:
                parsed_modules[mod_id] = MultiqcModule(mod_id, mod)
            else:
                # New sub-section
                parsed_modules[mod_id].update_init(c_id, mod)
            parsed_modules[mod_id].add_cc_section(c_id, mod)
            if mod["config"].get("plot_type") == "html":
                log.info(f"{c_id}: Found 1 sample (html)")
            elif mod["config"].get("plot_type") == "image":
                log.info(f"{c_id}: Found 1 sample (image)")
            else:
                log.info(f"{c_id}: Found {len(mod['data'])} samples ({mod['config'].get('plot_type')})")

    # Sort sections if we have a config option for order
    mod_order = getattr(config, "custom_content", {}).get("order", [])
    sorted_modules = [parsed_mod for parsed_mod in parsed_modules.values() if parsed_mod.anchor not in mod_order]
    sorted_modules.extend(
        [parsed_mod for mod_id in mod_order for parsed_mod in parsed_modules.values() if parsed_mod.anchor == mod_id]
    )

    # If we only have General Stats columns then there are no module outputs
    if len(sorted_modules) == 0:
        if len(cust_mod_by_id) == 1 and list(cust_mod_by_id.values())[0]["config"].get("plot_type") == "generalstats":
            sorted_modules = [bm]
        else:
            raise ModuleNoSamplesFound

    return sorted_modules


class MultiqcModule(BaseMultiqcModule):
    """Module class, used for each custom content type"""

    def __init__(self, c_id, mod):
        modname = c_id.replace("_", " ").title()
        mod_info = mod["config"].get("description")
        if "parent_name" in mod["config"]:
            modname = mod["config"]["parent_name"]
            mod_info = mod["config"].get("parent_description")
        elif "section_name" in mod["config"]:
            modname = mod["config"]["section_name"]
        if modname == "" or modname is None:
            modname = "Custom Content"

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name=modname,
            anchor=mod["config"].get("section_anchor", c_id),
            href=mod["config"].get("section_href"),
            info=mod_info,
            extra=mod["config"].get("extra"),
            # No DOI here.. // doi=
        )

        # Don't repeat the Custom Content name in the subtext
        if self.info or self.extra:
            self.intro = f"<p>{self.info}</p>{self.extra}"

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
        if self.info or self.extra:
            self.intro = f"<p>{self.info}</p>{self.extra}"

    def add_cc_section(self, c_id, mod):
        section_name = mod["config"].get("section_name", c_id.replace("_", " ").title())
        if section_name == "" or section_name is None:
            section_name = "Custom Content"

        section_description = mod["config"].get("description", "")

        pconfig = mod["config"].get("pconfig", {})
        if pconfig.get("id") is None:
            pconfig["id"] = f"{c_id}-plot"
        if pconfig.get("title") is None:
            pconfig["title"] = section_name

        plot = None
        content = None

        if not isinstance(mod["data"], list):
            mod["data"] = [mod["data"]]

        for i, ds in enumerate(mod["data"]):
            # Save the data if it's not a html string
            if not isinstance(ds, str):
                did = pconfig["id"]
                if i > 0:
                    did = f"{did}_{i}"
                self.write_data_file(ds, f"multiqc_{did}")
                pconfig["save_data_file"] = False

        # Try to coerce x-axis to numeric
        if mod["config"].get("plot_type") in ["linegraph", "scatter"]:
            try:
                mod["data"] = [{k: {float(x): v[x] for x in v} for k, v in ds.items()} for ds in mod["data"]]
            except ValueError:
                pass

        # Table
        if mod["config"].get("plot_type") == "table":
            pconfig["sort_rows"] = pconfig.get("sort_rows", pconfig.get("sortRows", False))
            headers = mod["config"].get("headers")
            plot = table.plot(mod["data"], headers, pconfig)

        # Bar plot
        elif mod["config"].get("plot_type") == "bargraph":
            mod["data"] = [{str(k): v for k, v in ds.items()} for ds in mod["data"]]
            plot = bargraph.plot(mod["data"], mod["config"].get("categories"), pconfig)

        # Line plot
        elif mod["config"].get("plot_type") == "linegraph":
            plot = linegraph.plot(mod["data"], pconfig)

        # Scatter plot
        elif mod["config"].get("plot_type") == "scatter":
            plot = scatter.plot(mod["data"], pconfig)

        # Heatmap
        elif mod["config"].get("plot_type") == "heatmap":
            plot = heatmap.plot(mod["data"], mod["config"].get("xcats"), mod["config"].get("ycats"), pconfig)

        # Violin plot
        elif mod["config"].get("plot_type") in ["violin", "beeswarm"]:
            plot = violin.plot(mod["data"], pconfig)

        # Raw HTML
        elif mod["config"].get("plot_type") == "html":
            if len(mod["data"]) > 1:
                log.warning(f"HTML plot type found with more than one dataset in {c_id}")
            content = mod["data"][0]

        # Raw image file as html
        elif mod["config"].get("plot_type") == "image":
            if len(mod["data"]) > 1:
                log.warning(f"Image plot type found with more than one dataset in {c_id}")
            content = mod["data"][0]

        # Not supplied
        elif mod["config"].get("plot_type") is None:
            log.warning(f"Plot type not found for content ID '{c_id}'")

        # Not recognised
        else:
            log.warning(
                "Error - custom content plot type '{}' not recognised for content ID {}".format(
                    mod["config"].get("plot_type"), c_id
                )
            )

        # Don't use exactly the same title / description text as the main module
        if section_name == self.name:
            section_name = None
        if self.info and section_description.strip(".") == self.info.strip("."):
            section_description = ""

        self.add_section(name=section_name, anchor=c_id, description=section_description, plot=plot, content=content)


def _find_file_header(f):
    # Collect commented out header lines
    hlines = []
    for line in f["f"].splitlines():
        if line.startswith("#"):
            hlines.append(line[1:])
    if len(hlines) == 0:
        return None
    hconfig = None
    try:
        hconfig = yaml.safe_load("\n".join(hlines))
        assert isinstance(hconfig, dict)
    except yaml.YAMLError as e:
        log.warning(f"Could not parse comment file header for MultiQC custom content: {f['fn']}")
        log.debug(e)
    except AssertionError:
        log.debug(f"Custom Content comment file header looked wrong: {hconfig}")
    else:
        return hconfig


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


def _parse_txt(f, conf):
    # Split the data into a list of lists by column
    sep = None
    if conf["file_format"] == "csv":
        sep = ","
    if conf["file_format"] == "tsv":
        sep = "\t"
    lines = f["f"].splitlines()
    d = []

    # Check for special case - HTML
    if conf.get("plot_type") == "html":
        for line in lines:
            if line and not line.startswith("#"):
                d.append(line)
        return "\n".join(d), conf

    # Not HTML, need to parse data
    ncols = None
    for line in lines:
        if line and not line.startswith("#"):
            sections = line.split(sep)
            d.append(sections)
            if ncols is None:
                ncols = len(sections)
            elif ncols != len(sections):
                log.warning(f"Inconsistent number of columns found in {f['fn']}! Skipping..")
                return None, conf

    # Convert values to floats if we can
    first_row_str = 0
    for i, line in enumerate(d):
        for j, v in enumerate(line):
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
            d[i][j] = v

    all_numeric = all([isinstance(line, float) for i in range(1, len(d)) for line in d[i][1:]])

    # General stat info files - expected to have at least 2 rows (first row always being the header)
    # and have at least 2 columns (first column always being sample name)
    if conf.get("plot_type") == "generalstats" and len(d) >= 2 and ncols >= 2:
        data = defaultdict(dict)
        for i, line in enumerate(d[1:], 1):
            for j, v in enumerate(line[1:], 1):
                data[line[0]][d[0][j]] = v
        return data, conf

    # Heatmap: Number of headers == number of lines
    if conf.get("plot_type") is None and first_row_str == len(lines) and all_numeric:
        conf["plot_type"] = "heatmap"
    if conf.get("plot_type") == "heatmap":
        conf["xcats"] = d[0][1:]
        conf["ycats"] = [s[0] for s in d[1:]]
        data = [s[1:] for s in d[1:]]
        return data, conf

    # Header row of strings, or configured as table
    if first_row_str == len(d[0]) or conf.get("plot_type") == "table":
        data = dict()
        for s in d[1:]:
            data[s[0]] = dict()
            for i, v in enumerate(s[1:]):
                cat = str(d[0][i + 1])
                data[s[0]][cat] = v
        # Bar graph or table - if numeric data, go for bar graph
        if conf.get("plot_type") is None:
            allfloats = True
            for r in d[1:]:
                for v in r[1:]:
                    allfloats = allfloats and isinstance(v, float)
            if allfloats:
                conf["plot_type"] = "bargraph"
            else:
                conf["plot_type"] = "table"
        # Set table col_1 header
        if conf.get("plot_type") == "table" and d[0][0].strip() != "":
            conf["pconfig"] = conf.get("pconfig", {})
            if not conf["pconfig"].get("col1_header"):
                conf["pconfig"]["col1_header"] = d[0][0].strip()
        # Return parsed data
        if conf.get("plot_type") == "bargraph" or conf.get("plot_type") == "table":
            return data, conf
        else:
            data = dict()  # reset

    # Scatter plot: First row is  str : num : num
    if (
        conf.get("plot_type") is None
        and len(d[0]) == 3
        and not isinstance(d[0][0], float)
        and isinstance(d[0][1], float)
        and isinstance(d[0][2], float)
    ):
        conf["plot_type"] = "scatter"

    if conf.get("plot_type") == "scatter":
        data = dict()
        for s in d:
            try:
                data[s[0]] = {"x": float(s[1]), "y": float(s[2])}
            except (IndexError, ValueError):
                pass
        return data, conf

    # Single sample line / bar graph - first row has two columns
    if len(d[0]) == 2:
        # Line graph - num : num
        if conf.get("plot_type") is None and isinstance(d[0][0], float) and isinstance(d[0][1], float):
            conf["plot_type"] = "linegraph"
        # Bar graph - str : num
        if conf.get("plot_type") is None and not isinstance(d[0][0], float) and isinstance(d[0][1], float):
            conf["plot_type"] = "bargraph"

        # Data structure is the same
        if conf.get("plot_type") == "linegraph" or conf.get("plot_type") == "bargraph":
            # Set section id based on directory if not known
            if conf.get("id") is None:
                conf["id"] = os.path.basename(f["root"])
            data = dict()
            for s in d:
                data[s[0]] = s[1]
            return {f["s_name"]: data}, conf

    # Multi-sample line graph: No header row, str : lots of num columns
    if conf.get("plot_type") is None and len(d[0]) > 4 and all_numeric:
        conf["plot_type"] = "linegraph"

    if conf.get("plot_type") == "linegraph":
        data = dict()
        # If the first row has no header, use it as axis labels
        x_labels = []
        if d[0][0].strip() == "":
            x_labels = d.pop(0)[1:]
        # Use 1..n range for x values
        for s in d:
            data[s[0]] = dict()
            for i, v in enumerate(s[1:]):
                try:
                    x_val = x_labels[i]
                    try:
                        x_val = float(x_val)
                    except ValueError:
                        pass
                except IndexError:
                    x_val = i + 1
                data[s[0]][x_val] = v
        return data, conf

    # Got to the end and haven't returned. It's a mystery, capn'!
    log.debug(
        f"Not able to figure out a plot type for '{f['fn']}' "
        + "plot type = {}, all numeric = {}, first row str = {}".format(
            conf.get("plot_type"), all_numeric, first_row_str
        )
    )
    return None, conf
