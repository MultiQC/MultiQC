import json
import logging
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, Union, List, Optional

from multiqc import report, config
from multiqc.base_module import BaseMultiqcModule
from multiqc.core.update_config import update_config, ClConfig
from multiqc.core.file_search import file_search
from multiqc.core.exec_modules import exec_modules
from multiqc.core.version_check import check_version
from multiqc.core.write_results import write_results
from multiqc.core.exceptions import RunError
from multiqc.plots.plotly.bar import BarPlot
from multiqc.plots.plotly.box import BoxPlot
from multiqc.plots.plotly.heatmap import HeatmapPlot
from multiqc.plots.plotly.line import LinePlot
from multiqc.plots.plotly.plot import PlotType, Plot
from multiqc.plots.plotly.scatter import ScatterPlot
from multiqc.plots.plotly.violin import ViolinPlot

# Set up logging
start_execution_time = time.time()
logger = logging.getLogger("multiqc")


def parse_logs(*analysis_dir, cfg: Optional[ClConfig] = None):
    """
    Parse files without generating a report. Useful to work with MultiQC interactively. Data can be accessed
    with other methods: `list_modules`, `show_plot`, `get_summarized_data`, etc.
    """
    update_config(*analysis_dir, cfg=cfg)

    check_version(parse_logs.__name__)

    report.reset_file_search()
    try:
        searched_modules = file_search()
        exec_modules(searched_modules, clean_up=False)
    except RunError as e:
        if e.message:
            logger.critical(e.message)


def parse_data_json(path: Union[str, Path]):
    """
    Try find multiqc_data.json in the given directory and load it into the report.
    """
    check_version(parse_data_json.__name__)

    json_path_found = False
    json_path: Path
    if path.endswith(".json"):
        json_path = path
        json_path_found = True
    else:
        json_path = Path(path) / "multiqc_data.json"
        if json_path.exists():
            json_path_found = True

    if not json_path_found:
        logger.error(f"multiqc_data.json not found in {path}")
        return

    # Loading from previous JSON
    logger.info(f"Loading data from {json_path}")
    try:
        with json_path.open("r") as f:
            data = json.load(f)

        for mod, sections in data["report_data_sources"].items():
            logger.info(f"Loaded module {mod}")
            for section, sources in sections.items():
                for sname, source in sources.items():
                    report.data_sources[mod][section][sname] = source
        for id, plot_dump in data["report_plot_data"].items():
            logger.info(f"Loaded plot {id}")
            report.plot_data[id] = plot_dump
    except (json.JSONDecodeError, KeyError) as e:
        logger.error(f"Error loading data from multiqc_data.json: {e}")


def list_data_sources() -> List[str]:
    """
    Return a list of the data sources that have been loaded.
    """
    file_list = []
    for mod, sections in report.data_sources.items():
        for section, sources in sections.items():
            for sname, source in sources.items():
                file_list.append(source)
    return file_list


def list_modules() -> List[str]:
    """
    Return a list of the modules that have been loaded in order.
    """
    return [m.name for m in report.modules_output]


def list_samples() -> List[str]:
    """
    Return a list of the samples that have been loaded.
    """
    samples = set()

    for mod, sections in report.data_sources.items():
        for section, sources in sections.items():
            for sname, source in sources.items():
                samples.add(sname)

    return sorted(samples)


def list_plots() -> Dict[str, List[Union[str, Dict[str, str]]]]:
    """
    Return a list of the plots that have been loaded for a given module,
    along with the number of datasets in each plot.
    """

    result = dict()
    for module in report.modules_output:
        result[module.name]: List[Union[str, Dict[str, str]]] = list()
        for section in module.sections:
            section_id = section["name"] or section["anchor"]
            if "plot_id" not in section:
                logger.warning(f"No plot_id in section {section_id} in module {module.name}")
                continue
            plot_id = section["plot_id"]
            if plot_id not in report.plot_by_id:
                raise ValueError(f'CRITICAL: Plot "{plot_id}" not found in report.plot_by_id')
            plot = report.plot_by_id[plot_id]
            if len(plot.datasets) == 1:
                result[module.name].append(section_id)
            if len(plot.datasets) > 1:
                result[module.name].append({section_id: [d.label for d in plot.datasets]})

    print(
        'List of available plots sections, by module. Use multiqc.show_plot("<module">, "<section>") to show plot. '
        + (
            "\nIf plot has several datasets, pass the dataset name: "
            'multiqc.show_plot("<module>", "<section>", "<dataset>")'
            if any(len(plot.datasets) > 1 for plot in report.plot_by_id.values())
            else ""
        )
    )
    return result


def show_plot(module: str, section: str, dataset: Optional[str] = None, **kwargs):
    """
    Show a plot in the notebook.
    """

    mod = next((m for m in report.modules_output if m.name == module or m.anchor == module), None)
    if not mod:
        raise ValueError(f'Module "{module}" is not found. Use multiqc.list_modules() to list available modules')

    sec = next((s for s in mod.sections if (s.name and s.name == section) or s.anchor == section), None)
    if not sec:
        raise ValueError(f'Section "{section}" is not found in module "{module}"')

    if sec.plot_id:
        plot = report.plot_by_id[sec.plot_id]
        ds_id = 0
        if dataset:
            for i, d in enumerate(plot.datasets):
                if d.label == dataset:
                    ds_id = i
                    break
        return plot.show(dataset_id=ds_id, **kwargs)
    elif sec.content:
        from IPython.core.display import HTML

        return HTML(sec.content)

    if dataset:
        raise ValueError(f'Plot section "{section}" with dataset "{dataset}" in module "{module}" not found')
    else:
        raise ValueError(f'Plot section "{section}" in module "{module}" not found')


def _load_plot(dump: Dict) -> Plot:
    """
    Load a plot and datasets from a JSON dump.
    """
    plot_type = PlotType(dump["plot_type"])
    if plot_type == PlotType.LINE:
        return LinePlot(**dump)
    elif plot_type == PlotType.BAR:
        return BarPlot(**dump)
    elif plot_type == PlotType.BOX:
        return BoxPlot(**dump)
    elif plot_type == PlotType.SCATTER:
        return ScatterPlot(**dump)
    elif plot_type == PlotType.HEATMAP:
        return HeatmapPlot(**dump)
    elif plot_type == PlotType.VIOLIN:
        return ViolinPlot(**dump)
    else:
        raise ValueError(f"Plot type {plot_type} is unknown or unsupported")


def get_general_stats_data(sample: Optional[str] = None) -> Dict:
    """
    Return parsed general stats data indexed by sample, then by data key. If sample is specified, return only data
    for that sample.
    """
    data = defaultdict(dict)
    for data_by_sample, header in zip(report.general_stats_data, report.general_stats_headers):
        for s, val_by_key in data_by_sample.items():
            if sample and s != sample:
                continue
            for key, val in val_by_key.items():
                if key in header:
                    key = f"{header[key].get('namespace', '')}.{key}"
                    data[s][key] = val
    if sample:
        if not data:
            return {}
        return data[sample]

    return data


def get_module_data(
    module: Optional[str] = None,
    sample: Optional[str] = None,
    key: Optional[str] = None,
) -> Dict:
    """
    Return parsed module data, indexed (optionally) by data key, then by sample. Module is either the module
    name, or the anchor.

    Takes data from report.saved_raw_data, which populated by self.write_data_file() calls in modules.
    This data is not necessarily normalized, e.g. numbers can be strings or numbers, depends on
    individual module behaviour.
    """
    if sample and sample not in list_samples():
        raise ValueError(f"Sample '{sample}' is not found. Use multiqc.list_samples() to list available samples")
    if module and module not in list_modules():
        raise ValueError(f"Module '{module}' is not found. Use multiqc.list_modules() to list available modules")

    data_by_module = {}
    for m in report.modules_output:
        if module and (m.name != module and m.anchor != module):
            continue

        module_data = m.saved_raw_data
        if sample:
            module_data = {k: v.get(sample, {}) for k, v in module_data.items()}
        if key:
            if module and key not in m.saved_raw_data:
                raise ValueError(f"Key '{key}' is not found in module '{module}'")
            module_data = module_data.get(key, {})
        elif len(module_data) == 1:  # only one key, flatten
            module_data = module_data[list(module_data.keys())[0]]

        data_by_module[m.name] = module_data

    if module:
        if not data_by_module:
            return {}
        return data_by_module[module]

    return data_by_module


def reset():
    """
    Reset the report to start fresh. Drops all previously parsed data.
    """
    config.reset()
    report.reset()


def add_custom_content_section(
    name,
    anchor,
    description="",
    content_before_plot="",
    plot: Optional[Union[Plot, str]] = None,
    content="",
    comment="",
    helptext="",
):
    """
    Add a custom content section to the report. This can be used to add a custom table or other content.
    """
    module = BaseMultiqcModule(
        name=name,
        anchor=anchor,
        info=description,
        comment=comment,
    )
    module.add_section(
        name=name,
        anchor=anchor,
        description=description,
        helptext=helptext,
        content_before_plot=content_before_plot,
        plot=plot,
        content=content,
        comment=comment,
    )
    report.modules_output.append(module)


def write_report(cfg: Optional[ClConfig] = None):
    """
    Write HTML and data files to disk. Useful to work with MultiQC interactively, after loading data with `load`.
    """
    update_config(cfg=cfg)

    check_version(write_report.__name__)

    if len(report.modules_output) == 0:
        logger.error("No analysis results found to make a report")
        return

    try:
        write_results()

    except RunError as e:
        if e.message:
            logger.critical(e.message)


def load_config(config_file: Union[str, Path]):
    """
    Load config on top of the current config from a MultiQC config file.
    """
    update_config()

    path = Path(config_file)
    if not path.exists():
        raise ValueError(f"Config file '{config_file}' not found")

    config.session_user_config_files.append(path.absolute())
    config.load_config_file(config_file)
