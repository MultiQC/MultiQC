""" Core MultiQC module to parse output from custom script output """


import base64
import json
import logging
import os
import re
from collections import OrderedDict, defaultdict

import yaml

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, beeswarm, heatmap, linegraph, scatter, table
from multiqc.utils import report

# Initialise the logger
log = logging.getLogger(__name__)


# Load YAML as an ordered dict
# From https://stackoverflow.com/a/21912744
def yaml_ordered_load(stream):
    class OrderedLoader(yaml.SafeLoader):
        pass

    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return OrderedDict(loader.construct_pairs(node))

    OrderedLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, construct_mapping)
    return yaml.load(stream, OrderedLoader)


def custom_module_classes():
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
    cust_mods = defaultdict(lambda: defaultdict(lambda: dict()))

    # Dictionary to hold search patterns - start with those defined in the config
    search_patterns = ["custom_content"]

    # First - find files using patterns described in the config
    config_data = getattr(config, "custom_data", {})
    mod_cust_config = {}
    for k, f in config_data.items():
        # Check that we have a dictionary
        if type(f) != dict:
            log.debug("config.custom_data row was not a dictionary: {}".format(k))
            continue
        c_id = f.get("id", k)

        # Data supplied in with config (e.g. from a multiqc_config.yaml file in working directory)
        if "data" in f:
            try:
                cust_mods[c_id]["data"].update(f["data"])
            except ValueError:
                # HTML plot type doesn't have a data sample-id key, so just take the whole chunk of data
                cust_mods[c_id]["data"] = f["data"]
            cust_mods[c_id]["config"].update({k: v for k, v in f.items() if k != "data"})
            cust_mods[c_id]["config"]["id"] = cust_mods[c_id]["config"].get("id", c_id)
            continue

        # Custom Content ID has search patterns in the config
        if c_id in report.files:
            cust_mods[c_id]["config"] = f
            cust_mods[c_id]["config"]["id"] = cust_mods[c_id]["config"].get("id", c_id)
            search_patterns.append(c_id)
            continue

        # Must just be configuration for a separate custom-content class
        mod_cust_config[c_id] = f

    # Now go through each of the file search patterns
    bm = BaseMultiqcModule()
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
                        parsed_data = yaml_ordered_load(f["f"])
                    except Exception as e:
                        log.warning("Error parsing YAML file '{}' (probably invalid YAML)".format(f["fn"]))
                        log.debug("YAML error: {}".format(e), exc_info=True)
                        break
                    parsed_data["id"] = parsed_data.get("id", f["s_name"])
                elif f_extension == ".json":
                    try:
                        # Use OrderedDict for objects so that column order is honoured
                        parsed_data = json.loads(f["f"], object_pairs_hook=OrderedDict)
                    except Exception as e:
                        log.warning("Error parsing JSON file '{}' (probably invalid JSON)".format(f["fn"]))
                        log.warning("JSON error: {}".format(e))
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
                    parsed_data.update(cust_mods.get(k, {}).get("config", {}))
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
                            cust_mods[c_id]["data"].update(parsed_data["data"])
                        else:
                            cust_mods[c_id]["data"] = parsed_data["data"]
                        cust_mods[c_id]["config"].update({j: k for j, k in parsed_data.items() if j != "data"})
                    else:
                        log.warning("No data found in {}".format(f["fn"]))

                # txt, csv, tsv etc
                else:
                    # Look for configuration details in the header
                    m_config = _find_file_header(f)
                    s_name = None
                    if m_config is not None:
                        c_id = m_config.get("id", k)
                        # Update the base config with anything parsed from the file
                        b_config = cust_mods.get(c_id, {}).get("config", {})
                        b_config.update(m_config)
                        # Now set the module config to the merged dict
                        m_config = dict(b_config)
                        s_name = m_config.get("sample_name")
                    else:
                        c_id = k
                        m_config = cust_mods.get(c_id, {}).get("config", {})

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
                            log.warning("Not able to parse custom data in {}".format(f["fn"]))
                        else:
                            # Did we get a new section id from the file?
                            if conf.get("id") is not None:
                                c_id = conf.get("id")
                            # heatmap - special data type
                            if type(parsed_data) == list:
                                cust_mods[c_id]["data"] = parsed_data
                            elif conf.get("plot_type") == "html":
                                cust_mods[c_id]["data"] = parsed_data
                            else:
                                cust_mods[c_id]["data"].update(parsed_data)
                            cust_mods[c_id]["config"].update(conf)
                    except (IndexError, AttributeError, TypeError):
                        log.error("Unexpected parsing error for {}".format(f["fn"]), exc_info=True)
                        raise  # testing
            except Exception as e:
                log.error("Uncaught exception raised for file '{}'".format(f["fn"]))
                log.exception(e)

        # Give log message if no files found for search pattern
        if num_sp_found_files == 0 and k != "custom_content":
            log.debug("No samples found: custom content ({})".format(k))

    # Filter to strip out ignored sample names
    for k in cust_mods:
        cust_mods[k]["data"] = bm.ignore_samples(cust_mods[k]["data"])

    # Remove any configs that have no data
    remove_cids = [k for k in cust_mods if len(cust_mods[k]["data"]) == 0]
    for k in remove_cids:
        del cust_mods[k]

    if len(cust_mods) == 0:
        raise UserWarning

    # Go through each data type
    parsed_modules = OrderedDict()
    for c_id, mod in cust_mods.items():
        # General Stats
        if mod["config"].get("plot_type") == "generalstats":
            gsheaders = mod["config"].get("pconfig")
            if gsheaders is None:
                headers = set()
                for d in mod["data"].values():
                    headers.update(d.keys())
                headers = list(headers)
                headers.sort()
                gsheaders = OrderedDict()
                for h in headers:
                    gsheaders[h] = dict()

            # Headers is a list of dicts
            if type(gsheaders) == list:
                gsheaders_dict = OrderedDict()
                for gsheader in gsheaders:
                    for col_id, col_data in gsheader.items():
                        gsheaders_dict[col_id] = col_data
                gsheaders = gsheaders_dict

            # Add namespace and description if not specified
            for m_id in gsheaders:
                if "namespace" not in gsheaders[m_id]:
                    gsheaders[m_id]["namespace"] = mod["config"].get("namespace", c_id)
            log.info("{}: Found {} General Statistics columns".format(c_id, len(mod["data"])))
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
                log.info("{}: Found 1 sample (html)".format(c_id))
            elif mod["config"].get("plot_type") == "image":
                log.info("{}: Found 1 sample (image)".format(c_id))
            else:
                log.info("{}: Found {} samples ({})".format(c_id, len(mod["data"]), mod["config"].get("plot_type")))

    # Sort sections if we have a config option for order
    mod_order = getattr(config, "custom_content", {}).get("order", [])
    sorted_modules = [parsed_mod for parsed_mod in parsed_modules.values() if parsed_mod.anchor not in mod_order]
    sorted_modules.extend(
        [parsed_mod for mod_id in mod_order for parsed_mod in parsed_modules.values() if parsed_mod.anchor == mod_id]
    )

    # If we only have General Stats columns then there are no module outputs
    if len(sorted_modules) == 0:
        if mod["config"].get("plot_type") == "generalstats":
            sorted_modules = [bm]
        else:
            raise UserWarning

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
            self.intro = "<p>{}</p>{}".format(self.info, self.extra)

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
            self.intro = "<p>{}</p>{}".format(self.info, self.extra)

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

        # Save the data if it's not a html string
        if not isinstance(mod["data"], str):
            self.write_data_file(mod["data"], "multiqc_{}".format(pconfig["id"]))
            pconfig["save_data_file"] = False

        # Try to cooerce x-axis to numeric
        if mod["config"].get("plot_type") in ["linegraph", "scatter"]:
            try:
                mod["data"] = {k: {float(x): v[x] for x in v} for k, v in mod["data"].items()}
            except ValueError:
                pass

        # Table
        if mod["config"].get("plot_type") == "table":
            pconfig["sortRows"] = pconfig.get("sortRows", False)
            headers = mod["config"].get("headers")
            plot = table.plot(mod["data"], headers, pconfig)

        # Bar plot
        elif mod["config"].get("plot_type") == "bargraph":
            mod["data"] = {k: v for k, v in sorted(mod["data"].items())}
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

        # Beeswarm plot
        elif mod["config"].get("plot_type") == "beeswarm":
            plot = beeswarm.plot(mod["data"], pconfig)

        # Raw HTML
        elif mod["config"].get("plot_type") == "html":
            content = mod["data"]

        # Raw image file as html
        elif mod["config"].get("plot_type") == "image":
            content = mod["data"]

        # Not supplied
        elif mod["config"].get("plot_type") == None:
            log.warning("Plot type not found for content ID '{}'".format(c_id))

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
    for l in f["f"].splitlines():
        if l.startswith("#"):
            hlines.append(l[1:])
    if len(hlines) == 0:
        return None
    hconfig = None
    try:
        hconfig = yaml.safe_load("\n".join(hlines))
        assert isinstance(hconfig, dict)
    except yaml.YAMLError as e:
        log.warning("Could not parse comment file header for MultiQC custom content: {}".format(f["fn"]))
        log.debug(e)
    except AssertionError:
        log.debug("Custom Content comment file header looked wrong: {}".format(hconfig))
    else:
        return hconfig


def _find_html_file_header(f):
    """Look for a HTML comment config at the start of a custom content HTML file"""
    if f["f"].lstrip().startswith("<!--"):
        match = re.search(r"^\<\!\-\-((?:.|\n|\r)*?)-->", f["f"].lstrip())
        if match:
            comment = match.group(1)
            try:
                return yaml_ordered_load(comment)
            except Exception as e:
                log.debug("Found Custom Content HTML comment, but couldn't load as YAML: {}".format(e), exc_info=True)
                log.debug("Comment:\n{}".format(comment))
    return {}


def _guess_file_format(f):
    """
    Tries to guess file format, first based on file extension (csv / tsv),
    then by looking for common column separators in the first 10 non-commented lines.
    Splits by tab / comma / space and counts resulting number of columns. Finds the most
    common column count, then comparsed how many lines had this number.
    eg. if tab, all 10 lines should have x columns when split by tab.
    Returns: csv | tsv | spaces   (spaces by default if all else fails)
    """
    filename, file_extension = os.path.splitext(f["fn"])
    tabs = []
    commas = []
    spaces = []
    j = 0
    for l in f["f"].splitlines():
        if not l.startswith("#"):
            j += 1
            tabs.append(len(l.split("\t")))
            commas.append(len(l.split(",")))
            spaces.append(len(l.split()))
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
        for l in lines:
            if l and not l.startswith("#"):
                d.append(l)
        return ("\n".join(d), conf)

    # Not HTML, need to parse data
    ncols = None
    for l in lines:
        if l and not l.startswith("#"):
            sections = l.split(sep)
            d.append(sections)
            if ncols is None:
                ncols = len(sections)
            elif ncols != len(sections):
                log.warning("Inconsistent number of columns found in {}! Skipping..".format(f["fn"]))
                return (None, conf)

    # Convert values to floats if we can
    first_row_str = 0
    for i, l in enumerate(d):
        for j, v in enumerate(l):
            if j != 0:  # we don't want to convert sample names to numbers
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

    all_numeric = all([type(l) == float for l in d[i][1:] for i in range(1, len(d))])

    # General stat info files - expected to have at least 2 rows (first row always being the header)
    # and have at least 2 columns (first column always being sample name)
    if conf.get("plot_type") == "generalstats" and len(d) >= 2 and ncols >= 2:
        data = defaultdict(dict)
        for i, l in enumerate(d[1:], 1):
            for j, v in enumerate(l[1:], 1):
                data[l[0]][d[0][j]] = v
        return (data, conf)

    # Heatmap: Number of headers == number of lines
    if conf.get("plot_type") is None and first_row_str == len(lines) and all_numeric:
        conf["plot_type"] = "heatmap"
    if conf.get("plot_type") == "heatmap":
        conf["xcats"] = d[0][1:]
        conf["ycats"] = [s[0] for s in d[1:]]
        data = [s[1:] for s in d[1:]]
        return (data, conf)

    # Header row of strings, or configured as table
    if first_row_str == len(d[0]) or conf.get("plot_type") == "table":
        data = OrderedDict()
        for s in d[1:]:
            data[s[0]] = OrderedDict()
            for i, v in enumerate(s[1:]):
                cat = str(d[0][i + 1])
                data[s[0]][cat] = v
        # Bar graph or table - if numeric data, go for bar graph
        if conf.get("plot_type") is None:
            allfloats = True
            for r in d[1:]:
                for v in r[1:]:
                    allfloats = allfloats and type(v) == float
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
            return (data, conf)
        else:
            data = OrderedDict()  # reset

    # Scatter plot: First row is  str : num : num
    if (
        conf.get("plot_type") is None
        and len(d[0]) == 3
        and type(d[0][0]) != float
        and type(d[0][1]) == float
        and type(d[0][2]) == float
    ):
        conf["plot_type"] = "scatter"

    if conf.get("plot_type") == "scatter":
        data = dict()
        for s in d:
            try:
                data[s[0]] = {"x": float(s[1]), "y": float(s[2])}
            except (IndexError, ValueError):
                pass
        return (data, conf)

    # Single sample line / bar graph - first row has two columns
    if len(d[0]) == 2:
        # Line graph - num : num
        if conf.get("plot_type") is None and type(d[0][0]) == float and type(d[0][1]) == float:
            conf["plot_type"] = "linegraph"
        # Bar graph - str : num
        if conf.get("plot_type") is None and type(d[0][0]) != float and type(d[0][1]) == float:
            conf["plot_type"] = "bargraph"

        # Data structure is the same
        if conf.get("plot_type") == "linegraph" or conf.get("plot_type") == "bargraph":
            # Set section id based on directory if not known
            if conf.get("id") is None:
                conf["id"] = os.path.basename(f["root"])
            data = OrderedDict()
            for s in d:
                data[s[0]] = s[1]
            return ({f["s_name"]: data}, conf)

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
        return (data, conf)

    # Got to the end and haven't returned. It's a mystery, capn'!
    log.debug(
        "Not able to figure out a plot type for '{}' ".format(f["fn"])
        + "plot type = {}, all numeric = {}, first row str = {}".format(
            conf.get("plot_type"), all_numeric, first_row_str
        )
    )
    return (None, conf)
