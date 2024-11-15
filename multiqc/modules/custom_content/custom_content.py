"""Core MultiQC module to parse output from custom script output"""

import base64
from io import BufferedReader
import json
import logging
import os
import re
from collections import defaultdict
from token import OP
from typing import Any, Dict, List, Optional, Set, Tuple, TypedDict, Union, cast

import yaml
from pydantic import BaseModel

from multiqc import Plot, config, report
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, box, heatmap, linegraph, scatter, table, violin
from multiqc.plots.plotly.bar import BarPlotConfig
from multiqc.plots.plotly.box import BoxPlotConfig
from multiqc.plots.plotly.heatmap import HeatmapConfig
from multiqc.plots.plotly.line import LinePlotConfig
from multiqc.plots.plotly.plot import PlotType
from multiqc.plots.plotly.scatter import ScatterConfig
from multiqc.plots.table_object import TableConfig
from multiqc.types import Anchor, LoadedFileDict, ModuleId, SectionId
from multiqc.validation import ModuleConfigValidationError

# Initialise the logger
log = logging.getLogger(__name__)


class CcDict(BaseModel):
    config: Dict[str, Any] = {}
    data: Union[Dict[str, Any], List[Dict[str, Any]], str] = {}


class ParsedDict(TypedDict, total=False):
    id: str
    plot_type: Optional[str]
    data: Union[str, Dict[str, Dict[str, Any]]]
    config: Dict[str, Any]
    section_name: Optional[str]


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
    ccdict_by_id: Dict[ModuleId, CcDict] = {}

    # Dictionary to hold search patterns - start with those defined in the config
    search_pattern_keys: List[ModuleId] = [ModuleId("custom_content")]

    # First - find files using patterns described in the config
    mod_cust_config = {}
    config_custom_data_id: ModuleId
    for config_custom_data_id, config_custom_data_item in config.custom_data.items():
        # Check that we have a dictionary
        if not isinstance(config_custom_data_item, dict):
            log.debug(f"config.custom_data row was not a dictionary: {config_custom_data_id}")
            continue
        cc_id: ModuleId = config_custom_data_item.get("id", config_custom_data_id)
        cc_mod_anchor: Anchor = config_custom_data_item.get("anchor", config_custom_data_id)

        # Data supplied in with config (e.g. from a multiqc_config.yaml file in working directory)
        if "data" in config_custom_data_item:
            ccdict: CcDict
            if cc_id in ccdict_by_id:
                ccdict = ccdict_by_id[cc_id]
                if isinstance(ccdict.data, dict):
                    assert isinstance(ccdict.data, dict)
                    ccdict.data.update(config_custom_data_item["data"])
                else:
                    # HTML plot type doesn't have a data sample-id key, so just take the whole chunk of data
                    ccdict.data = config_custom_data_item["data"]
            else:
                ccdict = CcDict(config={}, data=config_custom_data_item["data"])
                ccdict_by_id[cc_id] = ccdict

            assert isinstance(ccdict.config, dict)
            ccdict.config.update({k: v for k, v in config_custom_data_item.items() if k != "data"})
            if "id" not in ccdict.config:
                ccdict.config["id"] = cc_id
            if "anchor" not in ccdict.config:
                ccdict.config["anchor"] = cc_mod_anchor
            continue

        # Custom Content ID has search patterns in the config
        if report.files.get(cc_id):
            if "id" not in config_custom_data_item:
                config_custom_data_item["id"] = cc_id
            if "anchor" not in config_custom_data_item:
                config_custom_data_item["anchor"] = cc_mod_anchor
            if cc_id not in ccdict_by_id:
                ccdict_by_id[cc_id] = CcDict(config=config_custom_data_item)
            else:
                ccdict_by_id[cc_id].config = config_custom_data_item
            search_pattern_keys.append(cc_id)
            continue

        # Must just be configuration for a separate custom-content class
        mod_cust_config[cc_id] = config_custom_data_item

    bm: BaseMultiqcModule = BaseMultiqcModule(name="Custom content", anchor=Anchor("custom_content"))

    # Now go through each of the file search patterns
    for config_custom_data_id in search_pattern_keys:
        num_sp_found_files = 0
        for f in bm.find_log_files(config_custom_data_id):
            num_sp_found_files += 1

            f_extension = os.path.splitext(f["fn"])[1]

            # YAML and JSON files are the easiest
            parsed_dict: Optional[ParsedDict] = None
            if f_extension == ".yaml" or f_extension == ".yml":
                try:
                    parsed_dict = yaml.safe_load(f["f"])
                except Exception as e:
                    log.warning(f"Error parsing YAML file '{f['fn']}' (probably invalid YAML)")
                    log.debug(f"YAML error: {e}", exc_info=True)
                    break
                assert parsed_dict is not None
                parsed_dict["id"] = parsed_dict.get("id", f["s_name"])
            elif f_extension == ".json":
                try:
                    parsed_dict = json.loads(f["f"])
                except Exception as e:
                    log.warning(f"Error parsing JSON file '{f['fn']}' (probably invalid JSON)")
                    log.warning(f"JSON error: {e}")
                    break
                assert parsed_dict is not None
                parsed_dict["id"] = parsed_dict.get("id", f["s_name"])
            elif f_extension == ".png" or f_extension == ".jpeg" or f_extension == ".jpg":
                # image is an exception - will return a file handler
                image_string = base64.b64encode(cast(BufferedReader, f["f"]).read()).decode("utf-8")
                image_format = "png" if f_extension == ".png" else "jpg"
                img_html = '<div class="mqc-custom-content-image"><img src="data:image/{};base64,{}" /></div>'.format(
                    image_format, image_string
                )
                parsed_dict = {
                    "id": f["s_name"],
                    "plot_type": "image",
                    "section_name": f["s_name"].replace("_", " ").replace("-", " ").replace(".", " "),
                    "data": img_html,
                }
                # If the search pattern 'k' has an associated custom content section config, use it:
                if config_custom_data_id in ccdict_by_id:
                    parsed_dict.update(ccdict_by_id[config_custom_data_id].config)  # type: ignore
            elif f_extension == ".html":
                parsed_dict = {"id": f["s_name"], "plot_type": "html", "data": str(f["f"])}
                parsed_dict.update(_find_html_file_header(f))

            if parsed_dict is not None:
                parsed_item: Union[str, Dict, List, None] = parsed_dict.get("data", {})
                parsed_item_with_clean_sn: Dict[str, Dict[str, Any]] = {}
                if isinstance(parsed_item, dict):
                    # Run sample-name cleaning on the data keys
                    for sn, val_by_metric in parsed_item.items():
                        sn = bm.clean_s_name(sn, f)
                        if sn not in parsed_item_with_clean_sn:
                            parsed_item_with_clean_sn[sn] = val_by_metric
                        else:
                            parsed_item_with_clean_sn[sn].update(val_by_metric)
                    parsed_dict["data"] = parsed_item_with_clean_sn

                _c_id = parsed_dict.get("id", config_custom_data_id)
                parsed_item = parsed_dict.get("data", {})
                if parsed_item:
                    if _c_id not in ccdict_by_id:
                        ccdict_by_id[ModuleId(_c_id)] = CcDict()
                    _ccdict = ccdict_by_id[ModuleId(_c_id)]
                    if isinstance(parsed_item, dict) and isinstance(_ccdict.data, dict):
                        _ccdict.data.update(parsed_item)
                    else:
                        _ccdict.data = parsed_item
                    assert isinstance(_ccdict.config, dict)
                    _ccdict.config.update({j: k for j, k in parsed_dict.items() if j != "data"})
                else:
                    log.warning(f"No data found in {f['fn']}")

            # txt, csv, tsv etc
            else:
                # Look for configuration details in the header
                m_config: Optional[Dict[str, Any]]
                m_config, non_header_lines = _find_file_header(f)
                s_name = None
                c_id: ModuleId
                if m_config is not None:
                    c_id = m_config.get("id", config_custom_data_id)
                    # Update the base config with anything parsed from the file
                    b_config = ccdict_by_id.get(c_id, CcDict()).config
                    assert isinstance(b_config, dict)
                    b_config.update(m_config)
                    # Now set the module config to the merged dict
                    m_config = dict(b_config)
                    s_name = m_config.get("sample_name")
                else:
                    c_id = config_custom_data_id
                    m_config = (ccdict_by_id.get(c_id) or CcDict()).config

                # Guess sample name if not given
                if s_name is None:
                    s_name = f["s_name"]

                # Guess c_id if no information known
                if config_custom_data_id == "custom_content":
                    c_id = ModuleId(s_name)
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

                # Parse data and optionally guess plot type from data
                parsed_item, m_config, plot_type = _parse_txt(f, m_config, non_header_lines)
                m_config["plot_type"] = plot_type
                if parsed_item is None or len(parsed_item) == 0:
                    log.warning(f"Not able to parse custom data in {f['fn']}")
                else:
                    # Did we get a new section id from the file?
                    if m_config.get("id") is not None:
                        c_id = ModuleId(m_config["id"])
                    if c_id not in ccdict_by_id:
                        ccdict_by_id[c_id] = CcDict()
                    # heatmap - special data type
                    if isinstance(parsed_item, list):
                        ccdict_by_id[c_id].data = parsed_item
                    elif plot_type == PlotType.HTML:
                        ccdict_by_id[c_id].data = parsed_item
                    else:
                        assert isinstance(parsed_item, dict)
                        d = ccdict_by_id[c_id].data
                        assert isinstance(d, dict), (c_id, f["fn"], f["root"])
                        d.update(parsed_item)
                    assert isinstance(ccdict_by_id[c_id].config, dict)
                    ccdict_by_id[c_id].config.update(m_config)

        # Give log message if no files found for search pattern
        if num_sp_found_files == 0 and config_custom_data_id != "custom_content":
            log.debug(f"No samples found: custom content ({config_custom_data_id})")

    # Filter to strip out ignored sample names
    for config_custom_data_id in ccdict_by_id:
        ccdict = ccdict_by_id[config_custom_data_id]
        if isinstance(ccdict.data, dict):
            ccdict.data = bm.ignore_samples(ccdict.data)

    # Remove any configs that have no data
    remove_cids = [k for k in ccdict_by_id if len(ccdict_by_id[k].data) == 0]
    for config_custom_data_id in remove_cids:
        del ccdict_by_id[config_custom_data_id]

    if len(ccdict_by_id) == 0:
        raise ModuleNoSamplesFound

    # Go through each data type
    parsed_modules: Dict[ModuleId, MultiqcModule] = dict()
    mod_id: ModuleId
    for mod_id, ccdict in ccdict_by_id.items():
        # General Stats
        assert isinstance(ccdict.config, dict)
        plot_type = PlotType.from_str(ccdict.config.get("plot_type"))
        if plot_type == PlotType.GENERALSTATS:
            assert isinstance(ccdict.data, dict), ccdict.data
            gs_headers = ccdict.config.get("headers", ccdict.config.get("pconfig", {}))
            if gs_headers is None:
                headers_set: Set[str] = set()
                for _hd in ccdict.data.values():
                    headers_set.update(_hd.keys())
                headers = list(headers_set)
                headers.sort()
                gs_headers = dict()
                for h in headers:
                    gs_headers[h] = dict()

            # Headers is a list of dicts
            if isinstance(gs_headers, list):
                gs_headers_dict = dict()
                for gs_header in gs_headers:
                    for col_id, col_data in gs_header.items():
                        gs_headers_dict[col_id] = col_data
                gs_headers = gs_headers_dict

            # Add namespace and description if not specified
            for m_id in gs_headers:
                if "namespace" not in gs_headers[m_id]:
                    gs_headers[m_id]["namespace"] = ccdict.config.get("namespace", mod_id)
            log.info(f"{mod_id}: Found {len(ccdict.data)} General Statistics columns")
            bm.general_stats_addcols(ccdict.data, gs_headers)

        # Initialise this new module class and append to list
        else:
            # Is this file asking to be a sub-section under a parent section?
            mod_id = ccdict.config.get("parent_id", ccdict.config.get("id", mod_id))
            section_id: SectionId = ccdict.config.get("section_id", mod_id)

            mod_anchor: Optional[Anchor] = None
            if "parent_anchor" in ccdict.config:
                mod_anchor = ccdict.config["parent_anchor"]

            section_anchor: Optional[Anchor] = None
            if "section_anchor" in ccdict.config:
                section_anchor = ccdict.config["section_anchor"]
            elif "anchor" in ccdict.config:
                section_anchor = ccdict.config["anchor"]

            # Assuring anchors are unique, but prioritizing section over module
            if not section_anchor:
                section_anchor = Anchor(section_id)
            if not mod_anchor:
                mod_anchor = Anchor(mod_id)
            if section_anchor == mod_anchor:
                section_anchor = Anchor(f"{section_anchor}-section")

            # If we have any custom configuration from a MultiQC config file, update here
            # This is done earlier for tsv files too, but we do it here so that it overwrites what was in the file
            if mod_id in mod_cust_config:
                ccdict.config.update(mod_cust_config[mod_id])

            if mod_id in parsed_modules:
                # We've seen this module section before
                parsed_modules[mod_id].update_init(ccdict)
            else:
                parsed_modules[mod_id] = MultiqcModule(mod_id, anchor=mod_anchor, cc_dict=ccdict)

            parsed_modules[mod_id].add_cc_section(section_id, section_anchor=section_anchor, ccdict=ccdict)

            if plot_type == PlotType.HTML:
                log.info(f"{section_id}: Found 1 sample (html)")
            elif plot_type == PlotType.IMAGE:
                log.info(f"{section_id}: Found 1 sample (image)")
            else:
                log.info(f"{section_id}: Found {len(ccdict.data)} samples ({plot_type})")

    # Sort sections if we have a config option for order
    mod_order = getattr(config, "custom_content", {}).get("order", [])
    # after each element, also add a version with a "-module" next to it for back-compat with < 1.24
    mod_order = [x for x in mod_order for x in [x, re.sub("-module$", "", x), f"{x}-module"]]
    # remove duplicates, but keep order
    mod_order = list(dict.fromkeys(mod_order))
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
        for ccdict in ccdict_by_id.values():
            assert isinstance(ccdict.config, dict)
            cfgs.append(ccdict.config)
        if len(ccdict_by_id) >= 1 and all(
            PlotType.from_str(cfg.get("plot_type")) == PlotType.GENERALSTATS for cfg in cfgs
        ):
            sorted_modules = [bm]
        else:
            raise ModuleNoSamplesFound

    return sorted_modules


class MultiqcModule(BaseMultiqcModule):
    """Module class, used for each custom content type"""

    def __init__(self, id: ModuleId, anchor: Anchor, cc_dict: CcDict):
        modname: str = id.replace("_", " ").title()
        if modname == "":
            modname = "Custom Content"

        assert isinstance(cc_dict.config, dict)
        mod_info = cc_dict.config.get("description")
        if "parent_name" in cc_dict.config:
            assert isinstance(cc_dict.config["parent_name"], str)
            modname = cc_dict.config["parent_name"]
            mod_info = cc_dict.config.get("parent_description")
        elif "section_name" in cc_dict.config:
            assert isinstance(cc_dict.config["section_name"], str)
            modname = cc_dict.config["section_name"]

        super(MultiqcModule, self).__init__(
            name=modname,
            anchor=anchor,
            href=cc_dict.config.get("section_href"),
            info=mod_info,
            extra=cc_dict.config.get("extra"),
            doi=cc_dict.config.get("doi"),
        )
        self.id = id

        if "custom_content" in config.run_modules:
            # To allow file_search.include_or_exclude_modules() correctly filter these modules
            config.custom_content_modules.append(anchor)

    def update_init(self, ccdict: CcDict):
        """
        This function runs when we have already initialised a module
        and this is a subsequent file that will be another subsection.

        So most info should already be initialised properly.
        We check if anything is empty and if we have it set in this file
        we set it here.
        """

        if self.info is None or self.info == "":
            _cc_desc = ccdict.config.get("description")
            if _cc_desc is not None:
                self.info = str(_cc_desc)
        if self.extra is None or self.info == "":
            _cc_extra = ccdict.config.get("extra")
            if _cc_extra is not None:
                self.extra = str(_cc_extra)
        # This needs overwriting again as it has already run on init
        self.intro = self._get_intro()

    def add_cc_section(self, section_id: SectionId, section_anchor: Anchor, ccdict: CcDict):
        plot: Optional[Union[Plot[Any, Any], str]] = None
        content = None

        # Set section and plot id and anchor
        section_name: str = ccdict.config.get("section_name", section_id.replace("_", " ").title())
        if section_name == "":
            section_name = "Custom Content"

        pconfig = ccdict.config.get("pconfig", {})
        if pconfig.get("anchor") is None:
            if pconfig.get("id") is not None:
                pconfig["anchor"] = pconfig["id"]
                if pconfig["anchor"] == self.anchor or pconfig["anchor"] == section_anchor:
                    # making sure plot anchor is globally unique
                    pconfig["anchor"] += "-plot"
            else:
                pconfig["anchor"] = section_anchor + "-plot"  # making sure anchor is globally unique
        if pconfig.get("id") is None:
            pconfig["id"] = section_id
        if pconfig.get("title") is None:
            pconfig["title"] = section_name

        # But don't repeat header if it's the same title as the module title
        if section_name == self.name:
            section_name = ""

        plot_type: Optional[PlotType] = None
        _pt = ccdict.config.get("plot_type")
        if _pt is None:
            log.warning(f"'plot_type' value is not sepcified for content ID '{section_id}'")
        else:
            plot_type = PlotType.from_str(_pt)
            if plot_type is None:
                log.warning(f"Error: custom content plot type '{_pt}' not recognised for content ID {section_id}")

        plot_datasets: List[Any]  # to save after rendering

        # Heatmap
        if plot_type == PlotType.HEATMAP:
            plot = heatmap.plot(
                ccdict.data,  # type: ignore
                ccdict.config.get("xcats"),
                ccdict.config.get("ycats"),
                pconfig=HeatmapConfig(**pconfig),
            )
            plot_datasets = [ccdict.data]  # to save after rendering
        else:
            plot_datasets = [ccdict.data] if not isinstance(ccdict.data, list) else ccdict.data

            # Try to coerce x-axis to numeric
            if plot_type in [PlotType.LINE, PlotType.SCATTER]:
                try:
                    ccdict.data = [{k: {float(x): v[x] for x in v} for k, v in ds.items()} for ds in plot_datasets]
                except ValueError:
                    pass

            # Table
            if plot_type in [PlotType.TABLE, PlotType.VIOLIN]:
                headers = ccdict.config.get("headers")

                # handle some legacy fields for backwards compat
                sort_rows = pconfig.pop("sortRows", None)
                if sort_rows is not None:
                    pconfig["sort_rows"] = sort_rows
                no_violin = pconfig.pop("no_beeswarm", None)
                if no_violin is not None:
                    pconfig["no_violin"] = no_violin

                plot = (table if plot_type == PlotType.TABLE else violin).plot(
                    plot_datasets, headers=headers, pconfig=pconfig
                )

            # Bar plot
            elif plot_type == PlotType.BAR:
                ccdict.data = [{str(k): v for k, v in ds.items()} for ds in plot_datasets]
                plot = bargraph.plot(plot_datasets, ccdict.config.get("categories"), pconfig=BarPlotConfig(**pconfig))

            # Line plot
            elif plot_type == PlotType.LINE:
                plot = linegraph.plot(plot_datasets, pconfig=LinePlotConfig(**pconfig))  # type: ignore

            # Scatter plot
            elif plot_type == PlotType.SCATTER:
                plot = scatter.plot(ccdict.data, pconfig=ScatterConfig(**pconfig))  # type: ignore

            # Box plot
            elif plot_type == PlotType.BOX:
                plot = box.plot(plot_datasets, pconfig=BoxPlotConfig(**pconfig))  # type: ignore

            # Violin plot
            elif plot_type == PlotType.VIOLIN:
                plot = violin.plot(plot_datasets, pconfig=TableConfig(**pconfig))  # type: ignore

            # Raw HTML
            elif plot_type == PlotType.HTML:
                if isinstance(ccdict.data, list) and len(ccdict.data) > 1:
                    log.warning(f"HTML plot type found with more than one dataset in {section_id}")
                content = cast(str, ccdict.data)

            # Raw image file as html
            elif plot_type == PlotType.IMAGE:
                if isinstance(ccdict.data, list) and len(ccdict.data) > 1:
                    log.warning(f"Image plot type found with more than one dataset in {section_id}")
                content = cast(str, ccdict.data)

        if plot is not None:
            for i, ds in enumerate(plot_datasets):
                # Save the data if it's not a html string
                if not isinstance(plot, str):
                    did = plot.pconfig.id
                    if i > 0:
                        did = f"{did}_{i}"
                    self.write_data_file(ds, f"multiqc_{did}")
                    plot.pconfig.save_data_file = False

        if plot is not None or content:
            self.add_section(
                name=section_name,
                anchor=section_anchor,
                id=section_id,
                plot=plot,
                content=content or "",
            )


def _find_file_header(f: LoadedFileDict[str]) -> Tuple[Optional[Dict[str, Any]], List[str]]:
    # Collect commented out header lines
    hlines: List[str] = []
    other_lines: List[str] = []
    line: str
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
        raise ModuleConfigValidationError(
            f"Could not parse comment file header for MultiQC custom content: {f['fn']}. "
            + "Note that everything behind a comment character '#' is expected to in YAML format. "
            + "To provide column names in a TSV or CSV file, put them as the first raw without fencing it with a '#'.",
            "custom_content",
        )
    if not isinstance(hconfig, dict):
        raise ModuleConfigValidationError(
            "Custom Content comment file header looked wrong. It's expected to "
            + f"be parsed to a dict, got {type(hconfig)}: {hconfig}",
            "custom_content",
        )
    return cast(Dict[str, Any], hconfig), other_lines


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


def _parse_txt(
    f, conf: Dict, non_header_lines: List[str]
) -> Tuple[Union[str, Dict, List, None], Dict, Optional[PlotType]]:
    """
    Parse data and optionally guess plot type
    """

    # Split the data into a list of lists by column
    sep = None
    if conf.get("file_format") == "csv":
        sep = ","
    if conf.get("file_format") == "tsv":
        sep = "\t"

    plot_type: Optional[PlotType] = None
    _pt: Optional[str] = conf.get("plot_type")
    if _pt is not None:
        plot_type = PlotType.from_str(_pt)

    # Check for special case - HTML
    if plot_type == PlotType.HTML:
        lines: List[str] = []
        for line in non_header_lines:
            if line:
                lines.append(line)
        return "\n".join(lines), conf, plot_type

    # Not HTML, need to parse data
    matrix: List[List[Union[str, float, int]]] = []
    ncols = None
    for line in non_header_lines:
        if line.rstrip():
            sections = cast(List[Any], line.rstrip().split(sep))
            matrix.append(sections)
            if ncols is None:
                ncols = len(sections)
            elif ncols != len(sections):
                log.warning(f"Inconsistent number of columns found in {f['fn']}! Skipping..")
                return None, conf, plot_type

    # Convert values to floats if we can
    strings_in_first_row = 0
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
                    strings_in_first_row += 1
            matrix[i][j] = v

    all_numeric = all(isinstance(v, (float, int)) for s in matrix[1:] for v in s[1:])

    # General stat info files - expected to have at least 2 rows (first row always being the header)
    # and have at least 2 columns (first column always being sample name)
    if plot_type == PlotType.GENERALSTATS and len(matrix) >= 2 and ncols and ncols >= 2:
        data_ddict: Dict[str, Dict] = defaultdict(dict)
        for i, section in enumerate(matrix[1:], 1):
            for j, v in enumerate(section[1:], 1):
                s_name = str(section[0])
                data_ddict[s_name][matrix[0][j]] = v
        return data_ddict, conf, plot_type

    # Heatmap: number of headers == number of lines
    if plot_type is None and strings_in_first_row == len(non_header_lines) and all_numeric:
        plot_type = PlotType.HEATMAP
    if plot_type == PlotType.HEATMAP:
        conf["xcats"] = matrix[0][1:]
        conf["ycats"] = [s[0] for s in matrix[1:]]
        data_list: List = [s[1:] for s in matrix[1:]]
        return data_list, conf, plot_type

    # Header row of strings box plot
    if strings_in_first_row == len(matrix[0]) and plot_type == PlotType.BOX:
        box_ddict: Dict[str, List[float]] = {str(v): [] for v in matrix[0]}
        for sidx, s in enumerate(matrix[0]):
            for row in matrix[1:]:
                try:
                    val = float(row[sidx])
                except ValueError:
                    pass
                else:
                    box_ddict[str(s)].append(val)
        return box_ddict, conf, plot_type

    # Header row of strings, or configured as table
    if strings_in_first_row == len(matrix[0]) or plot_type in [PlotType.TABLE, PlotType.VIOLIN]:
        data_ddict = dict()
        for row in matrix[1:]:
            s_name = str(row[0])
            data_ddict[s_name] = dict()
            for i, v in enumerate(row[1:]):
                cat = str(matrix[0][i + 1])
                data_ddict[s_name][cat] = v
        # Bar graph or table - if numeric data, go for bar graph
        if plot_type is None:
            allfloats = True
            for r in matrix[1:]:
                for v in r[1:]:
                    allfloats = allfloats and isinstance(v, float)
            if allfloats:
                plot_type = PlotType.BAR
            else:
                plot_type = PlotType.TABLE
        # Set table col_1 header
        col_name = str(matrix[0][0])
        if plot_type in [PlotType.TABLE, PlotType.VIOLIN] and col_name.strip() != "":
            conf["pconfig"] = conf.get("pconfig", {})
            if not conf["pconfig"].get("col1_header"):
                conf["pconfig"]["col1_header"] = col_name.strip()
        # Return parsed data
        if plot_type in [PlotType.BAR, PlotType.TABLE, PlotType.VIOLIN]:
            return data_ddict, conf, plot_type
        else:
            data_ddict = dict()  # reset

    # Scatter plot: first row is str : num : num
    if (
        plot_type is None
        and len(matrix[0]) == 3
        and not isinstance(matrix[0][0], float)
        and isinstance(matrix[0][1], float)
        and isinstance(matrix[0][2], float)
    ):
        plot_type = PlotType.SCATTER

    if plot_type == PlotType.SCATTER:
        dicts: Dict[str, Dict[str, float]] = dict()
        for row in matrix:
            try:
                dicts[str(row[0])] = {"x": float(row[1]), "y": float(row[2])}
            except (IndexError, ValueError):
                pass
        return dicts, conf, plot_type

    # Single-sample box/bar/line plot
    if len(matrix[0]) == 2 or len(matrix[0]) == 1:
        # Set section id based on directory if not known
        if conf.get("id") is None:
            conf["id"] = os.path.basename(f["root"])

        # Single-sample box plot
        if len(matrix[0]) == 1:
            if plot_type is None and isinstance(matrix[0][0], float):
                plot_type = PlotType.BOX
            if plot_type == PlotType.BOX:
                return {f["s_name"]: [float(r[0]) for r in matrix]}, conf, plot_type

        # Single sample line/bar plot - first row has two columns
        if len(matrix[0]) == 2:
            # Line graph - num : num
            if plot_type is None and isinstance(matrix[0][0], float) and isinstance(matrix[0][1], float):
                plot_type = PlotType.LINE
            # Bar graph - str : num
            if plot_type is None and not isinstance(matrix[0][0], float) and isinstance(matrix[0][1], float):
                plot_type = PlotType.BAR

            # Data structure is the same
            if plot_type in [PlotType.LINE, PlotType.BAR]:
                data_dict: Dict = dict()
                for row in matrix:
                    data_dict[row[0]] = row[1]
                return {f["s_name"]: data_dict}, conf, plot_type

    # Multi-sample line graph: No header row, lots of num columns
    if plot_type is None and len(matrix[0]) > 4 and all_numeric:
        plot_type = PlotType.LINE

    if plot_type == PlotType.LINE:
        data_ddict = dict()
        # If the first row has no header, use it as axis labels
        x_labels = []
        s_name = str(matrix[0][0])
        if s_name.strip() == "":
            x_labels = matrix.pop(0)[1:]
        # Use 1..n range for x values
        for row in matrix:
            name = str(row[0])
            data_ddict[name] = dict()
            for i, v in enumerate(row[1:]):
                try:
                    x_val = x_labels[i]
                    try:
                        x_val = float(x_val)
                    except ValueError:
                        pass
                except IndexError:
                    x_val = i + 1
                data_ddict[name][x_val] = v
        return data_ddict, conf, plot_type

    # Got to the end and haven't returned. It's a mystery, capn'!
    log.debug(
        f"Not able to figure out a plot type for '{f['fn']}' "
        + f"plot type = {plot_type}, all numeric = {all_numeric}, first row str = {strings_in_first_row}"
    )
    return None, conf, plot_type
