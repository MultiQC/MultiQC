"""
Special case MultiQC module to load multiqc.parquet
It allows rerunning MultiQC when original data is gone, as well as extend
existing reports with new data.
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import packaging
import pandas as pd

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


def load_plot_input(plot_df: pd.DataFrame) -> Tuple[NormalizedPlotInputData, Union[Plot, str, None]]:
    plot_type = PlotType.from_str(plot_df.plot_type.iloc[0])
    pconfig_str = plot_df._pconfig.iloc[0]
    if isinstance(pconfig_str, str):
        pconfig = json.loads(pconfig_str)
    anchor = Anchor(str(plot_df.anchor.iloc[0]))

    plot_input: NormalizedPlotInputData
    plot: Union[Plot, str, None]

    if plot_type == PlotType.LINE:
        plot_input = LinePlotNormalizedInputData.from_df(plot_df, pconfig, anchor)
        plot = LinePlot.from_inputs(plot_input)
    elif plot_type == PlotType.BAR:
        plot_input = BarPlotInputData.from_df(plot_df, pconfig, anchor)
        plot = BarPlot.from_inputs(plot_input)
    elif plot_type == PlotType.BOX:
        plot_input = BoxPlotInputData.from_df(plot_df, pconfig, anchor)
        plot = BoxPlot.from_inputs(plot_input)
    elif plot_type == PlotType.HEATMAP:
        plot_input = HeatmapNormalizedInputData.from_df(plot_df, pconfig, anchor)
        plot = HeatmapPlot.from_inputs(plot_input)
    elif plot_type == PlotType.VIOLIN or plot_type == PlotType.TABLE:
        plot_input = ViolinPlotInputData.from_df(plot_df, pconfig, anchor)
        plot = ViolinPlot.from_inputs(plot_input)
    elif plot_type == PlotType.SCATTER:
        plot_input = ScatterNormalizedInputData.from_df(plot_df, pconfig, anchor)
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
            # Extract metadata from parquet
            metadata = plot_data_store.get_report_metadata(path)
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

                    # Get intro and comment
                    intro = mod_dict.pop("intro", "")
                    comment = mod_dict.pop("comment", "")

                    # Create module
                    mod = BaseMultiqcModule(name=name, anchor=Anchor(anchor), info=info)
                    mod.sections = sections
                    mod.versions = versions
                    mod.intro = intro
                    mod.comment = comment

                    log.info(f"Loading module {mod.name} from parquet")
                    report.modules.append(mod)

            # Load data sources
            if "data_sources" in metadata:
                for mod_id, source_dict in metadata["data_sources"].items():
                    for section, sources in source_dict.items():
                        for sname, source in sources.items():
                            report.data_sources[mod_id][section][sname] = source

            # Set creation date
            if "creation_date" in metadata:
                try:
                    report.creation_date = datetime.strptime(metadata["creation_date"], "%Y-%m-%d, %H:%M %Z")
                except ValueError:
                    try:
                        report.creation_date = datetime.strptime(metadata["creation_date"], "%Y-%m-%d, %H:%M")
                    except ValueError:
                        log.error(f"Could not parse creation date: {metadata['creation_date']}")

            # Load config values
            if "config" in metadata:
                log.debug("Loading config values from parquet")
                # TOOD: Update loaded configs

            df = pd.read_parquet(path)
            for anchor, plot_df in df.groupby("anchor"):
                if "plot_type" in plot_df.columns and plot_df.plot_type.iloc[0]:
                    anchor = Anchor(str(anchor))
                    plot_input, plot = load_plot_input(plot_df)
                    report.plot_input_data[anchor] = plot_input
                    if plot is not None:
                        report.plot_by_id[anchor] = plot

        except Exception as e:
            log.error(f"Error loading data from parquet file: {e}")
            if config.strict:
                raise e
