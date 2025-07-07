"""
Special case MultiQC module to load multiqc.parquet
It allows rerunning MultiQC when original data is gone, as well as extend
existing reports with new data.
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import packaging.version
import polars as pl

from multiqc import config, report
from multiqc.base_module import BaseMultiqcModule, Section
from multiqc.core import plot_data_store
from multiqc.plots.bargraph import BarPlot, BarPlotInputData
from multiqc.plots.box import BoxPlot, BoxPlotInputData
from multiqc.plots.heatmap import HeatmapNormalizedInputData, HeatmapPlot
from multiqc.plots.linegraph import LinePlot, LinePlotNormalizedInputData
from multiqc.plots.plot import NormalizedPlotInputData, Plot
from multiqc.plots.scatter import ScatterNormalizedInputData, ScatterPlot
from multiqc.plots.violin import ViolinPlot, ViolinPlotInputData
from multiqc.types import Anchor, PlotType

log = logging.getLogger(__name__)


def load_plot_input(plot_input_data_dict: Dict) -> Tuple[NormalizedPlotInputData, Union[Plot, str, None]]:
    # Process the JSON data to replace NaN markers with proper NaN values
    def replace_nan_markers(obj: Any) -> Any:
        import math

        if isinstance(obj, str) and obj == "__NAN__MARKER__":
            return math.nan
        elif isinstance(obj, dict):
            return {k: replace_nan_markers(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [replace_nan_markers(i) for i in obj]
        return obj

    # Ensure we're working with a dictionary and apply the replacement to each value
    if not isinstance(plot_input_data_dict, dict):
        raise TypeError(f"Expected a dictionary, got {type(plot_input_data_dict)}")

    processed_dict = {k: replace_nan_markers(v) for k, v in plot_input_data_dict.items()}
    plot_type = PlotType.from_str(processed_dict["plot_type"])

    plot_input: NormalizedPlotInputData
    plot: Union[Plot, str, None]

    if plot_type == PlotType.LINE:
        plot_input = LinePlotNormalizedInputData(**processed_dict)
        plot = LinePlot.from_inputs(plot_input)
    elif plot_type == PlotType.BAR:
        plot_input = BarPlotInputData(**processed_dict)
        plot = BarPlot.from_inputs(plot_input)
    elif plot_type == PlotType.BOX:
        plot_input = BoxPlotInputData(**processed_dict)
        plot = BoxPlot.from_inputs(plot_input)
    elif plot_type == PlotType.HEATMAP:
        plot_input = HeatmapNormalizedInputData(**processed_dict)
        plot = HeatmapPlot.from_inputs(plot_input)
    elif plot_type == PlotType.VIOLIN or plot_type == PlotType.TABLE:
        plot_input = ViolinPlotInputData(**processed_dict)
        plot = ViolinPlot.from_inputs(plot_input)
    elif plot_type == PlotType.SCATTER:
        plot_input = ScatterNormalizedInputData(**processed_dict)
        plot = ScatterPlot.from_inputs(plot_input)
    else:
        raise ValueError(f"Unknown plot type: {plot_type}")

    return plot_input, plot


class LoadMultiqcData(BaseMultiqcModule):
    def __init__(self):
        super(LoadMultiqcData, self).__init__(
            name="MultiQC Data",
            anchor=Anchor("multiqc_data"),
            info="loads multiqc data",
        )

        # First, try to find parquet file
        parquet_files = self.find_log_files("multiqc_data")
        if parquet_files:
            for f in parquet_files:
                self.load_parquet_file(Path(f["root"]) / f["fn"])

    def load_parquet_file(self, path: Union[str, Path]):
        """
        Load a multiqc.parquet file containing all report data.
        """
        path = Path(path)
        assert path.suffix == ".parquet"
        log.info(f"Loading report data from parquet file: {path}")

        try:
            # Read the entire parquet file
            df = pl.read_parquet(path)

            # Extract metadata from parquet
            metadata = plot_data_store.get_report_metadata(df)
            if metadata is None:
                log.error(f"Failed to extract metadata from parquet file: {path}")
                return

            # Load modules
            if "modules" in metadata:
                for mod_dict in metadata["modules"]:
                    sections = [Section(**section) for section in mod_dict.pop("sections")]

                    # Convert versions to expected format
                    versions: Dict[str, List[Tuple[Optional[packaging.version.Version], str]]] = {}
                    if "versions" in mod_dict:
                        versions_data = mod_dict.pop("versions")
                        versions = {
                            name: [(None, version) for version in versions] for name, versions in versions_data.items()
                        }

                    # Extract other module data
                    anchor = mod_dict.pop("anchor")
                    name = mod_dict.pop("name")
                    info = mod_dict.pop("info", "")
                    intro = mod_dict.pop("intro", "")
                    comment = mod_dict.pop("comment", "")

                    # Special handling for Software Versions modules - skip them to avoid duplicates
                    if anchor == "multiqc_software_versions" and name == "Software Versions":
                        # Extract software versions from the HTML content
                        if sections:
                            log.info("Extracting software versions from HTML content in parquet file")
                            import re

                            for section in sections:
                                if section.content:
                                    # Parse the HTML table to extract software versions
                                    # Look for table rows with pattern: <td>Software</td><td><samp>Version</samp></td>
                                    row_pattern = r"<tr><td>([^<]+)</td><td><samp>([^<]+)</samp></td></tr>"
                                    matches = re.findall(row_pattern, section.content)

                                    for software_name, version in matches:
                                        # Add to global software_versions using software name as group
                                        report.software_versions[software_name][software_name].append(version)
                                        log.debug(f"Added software version: {software_name} = {version}")

                        log.info(
                            "Skipping Software Versions module from parquet file - extracted data to global software_versions"
                        )
                        continue  # Skip adding this module to report.modules

                    # Create module
                    mod = BaseMultiqcModule(name=name, anchor=Anchor(anchor), info=info)
                    mod.sections = sections
                    mod.versions = versions
                    mod.intro = intro
                    mod.comment = comment

                    log.info(f"Loading module {mod.name} from parquet")
                    report.modules.append(mod)

            # Load global software versions data
            if "software_versions" in metadata:
                software_versions_data = metadata["software_versions"]
                log.info("Loading global software versions data from parquet file")
                for group_name, group_versions in software_versions_data.items():
                    for software_name, versions_list in group_versions.items():
                        for version_str in versions_list:
                            report.software_versions[group_name][software_name].append(version_str)

            # Load data sources
            if "data_sources" in metadata:
                for mod_id, source_dict in metadata["data_sources"].items():
                    for section, sources in source_dict.items():
                        for sname, source in sources.items():
                            report.data_sources[mod_id][section][sname] = source

            # Set creation date
            if "creation_date" in metadata:
                try:
                    # Convert the datetime to a Python datetime object
                    if isinstance(metadata["creation_date"], datetime):
                        report.creation_date = metadata["creation_date"]
                    else:
                        creation_date_str = str(metadata["creation_date"])
                        # Use standard Python datetime parsing
                        report.creation_date = datetime.fromisoformat(creation_date_str.replace("Z", "+00:00"))
                except ValueError as e:
                    log.error(f"Could not parse creation date: {metadata['creation_date']}, error: {e}")

            if "config" in metadata:
                pass  # We do not load config, but keep the current one

            # Load plot input data from plot_input rows in the dataframe
            if "type" in df.columns and "plot_input_data" in df.columns:
                plot_input_rows = df.filter(pl.col("type") == "plot_input")
                for row in plot_input_rows.iter_rows(named=True):
                    anchor = Anchor(str(row["anchor"]))
                    plot_input_data = row["plot_input_data"]
                    plot_input_data_dict = json.loads(plot_input_data)
                    try:
                        plot_input, plot = load_plot_input(plot_input_data_dict)
                        report.plot_input_data[anchor] = plot_input
                        if plot is not None:
                            report.plot_by_id[anchor] = plot
                    except Exception as e:
                        log.error(f"Error loading plot input data {anchor}: {e}")
                        if config.strict:
                            raise e

        except Exception as e:
            log.error(f"Error loading data from parquet file: {e}")
            if config.strict:
                raise e
