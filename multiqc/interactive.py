import json
import logging
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, Union, List, Optional

from multiqc import core
from multiqc.core import RunError
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots.plotly.bar import BarPlot
from multiqc.plots.plotly.box import BoxPlot
from multiqc.plots.plotly.heatmap import HeatmapPlot
from multiqc.plots.plotly.line import LinePlot
from multiqc.plots.plotly.plot import go, PlotType, Plot
from multiqc.plots.plotly.scatter import ScatterPlot
from multiqc.plots.plotly.violin import ViolinPlot
from multiqc.utils import report, config
from multiqc.utils.copy_function_signature import copy_callable_signature

# Set up logging
start_execution_time = time.time()
logger = logging.getLogger("multiqc")


def load_config(config_file: str):
    """
    Load config from a MultiQC config file.
    """
    config.user_config_files.append(Path(config_file).absolute())
    config.load_config(config_file)


@copy_callable_signature(core.init_config)
def parse_logs(analysis_dir: Union[str, List[str]], **kwargs):
    """
    Parse files without generating a report. Useful to work with MultiQC interactively. Data can be accessed
    with other methods: `list_modules`, `show_plot`, `get_summarized_data`, etc.
    """
    core.init_config(analysis_dir, **kwargs)

    # We want to keep report session, so we can interactively append modules to the report
    # if not report.initialized:
    #     report.reset()

    report.reset_file_search()

    try:
        searched_modules = core.file_search()

        core.exec_modules(searched_modules, clean_up=False)

    except RunError as e:
        if e.message:
            logger.critical(e.message)


def parse_data_json(analysis_dir: Union[str, List[str]]):
    """
    Try find multiqc_data.json in the given directory and load it into the report.
    """

    json_path_found = False
    json_path = None
    if not isinstance(analysis_dir, list):
        json_path = Path("multiqc_data.json")
        if json_path.exists():
            json_path_found = True
        else:
            json_path = Path(analysis_dir) / "multiqc_data.json"
            if json_path.exists():
                json_path_found = True

    if not json_path_found:
        logger.error("multiqc_data.json not found in the given directory")

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


def list_plots() -> List[str]:
    """
    Return a list of the plots that have been loaded for a given module,
    along with the number of datasets in each plot.
    """
    return list(
        f"{id} ({len(plot['datasets'])} datasets)" if len(plot["datasets"]) > 1 else id
        for id, plot in report.plot_data.items()
    )


def _load_plot(dump: Dict) -> Plot:
    """
    Load a plot and datasets from a JSON dump.
    """
    plot_type = PlotType(dump["plot_type"])
    dump["layout"] = go.Layout(dump["layout"])
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


def show_plot(plot_id: str, dataset_id=0, **kwargs) -> go.Figure:
    """
    Show a plot in the notebook.
    """
    dump = report.plot_data[plot_id]
    model = _load_plot(dump)
    return model.show(dataset_id=dataset_id, **kwargs)


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
    # report.plot_data[self.id] = self.model_dump(warnings=False)
    # report.plot_data[self.id]["layout"] = self.layout.to_plotly_json()
    report.modules_output.append(module)


@copy_callable_signature(core.init_config)
def write_report(**kwargs):
    """
    Write HTML and data files to disk. Useful to work with MultiQC interactively, after loading data with `load`.
    """
    core.init_config(**kwargs)

    if len(report.modules_output) == 0:
        logger.error("No analysis results found to make a report")
        return

    try:
        core.write_results()

    except RunError as e:
        if e.message:
            logger.critical(e.message)
