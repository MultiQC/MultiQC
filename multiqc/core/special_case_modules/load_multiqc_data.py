"""Special case MultiQC module to load multiqc_data.json
It allows rerunning MultiQC when original data is gone, as well as extend
existing reports with new data.
"""

import json
import logging
import os
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, cast

import packaging.version
import pandas as pd

from multiqc import report
from multiqc.base_module import BaseMultiqcModule, Section
from multiqc.core import tmp_dir
from multiqc.plots.bargraph import BarPlot
from multiqc.plots.box import BoxPlot
from multiqc.plots.heatmap import HeatmapPlot
from multiqc.plots.linegraph import LinePlot
from multiqc.plots.plot import NormalizedPlotInputData, Plot
from multiqc.plots.scatter import ScatterPlot
from multiqc.plots.table_object import ColumnKey, InputRow, SampleGroup, SampleName
from multiqc.plots.violin import ViolinPlot
from multiqc.types import Anchor, PlotType

log = logging.getLogger(__name__)


# def load_plot(plot_dump: Dict) -> Plot:
#     if plot_dump["plot_type"] == PlotType.LINE.value:
#         return LinePlot(**plot_dump)
#     elif plot_dump["plot_type"] == PlotType.BAR.value:
#         # We missing category with a "nan". however, json doesn't allow nans,
#         # so it is replaced with null in dump. so we need to put nans back
#         for ds in plot_dump["datasets"]:
#             for cat in ds["cats"]:
#                 cat["data"] = ["nan" if x is None else x for x in cat["data"]]
#                 cat["data_pct"] = ["nan" if x is None else x for x in cat["data_pct"]]
#         return BarPlot(**plot_dump)
#     elif plot_dump["plot_type"] == PlotType.BOX.value:
#         return BoxPlot(**plot_dump)
#     elif plot_dump["plot_type"] == PlotType.HEATMAP.value:
#         return HeatmapPlot(**plot_dump)
#     elif plot_dump["plot_type"] == PlotType.VIOLIN.value:
#         return ViolinPlot(**plot_dump)
#     elif plot_dump["plot_type"] == PlotType.SCATTER.value:
#         return ScatterPlot(**plot_dump)
#     else:
#         raise ValueError(f"Unknown plot type: {plot_dump['plot_type']}")


class LoadMultiqcData(BaseMultiqcModule):
    def __init__(self):
        super(LoadMultiqcData, self).__init__(
            name="MultiQC Data",
            anchor=Anchor("multiqc_data"),
            info="loads multiqc_data.json",
        )

        for f in self.find_log_files("multiqc_data"):
            self.load_data_json(Path(f["root"]) / f["fn"])

    def load_data_json(self, path: Union[str, Path]):
        """
        Try find multiqc_data.json in the given directory, and load it into the report.
        """
        path = Path(path)
        assert path.suffix == ".json"
        log.info(f"Loading previous run from {path}")
        try:
            with path.open("r") as f:
                data = json.load(f)

            # Load module instances (doesn't include data)
            for mod_dict in data["report_modules"]:
                sections = [Section(**section) for section in mod_dict.pop("sections")]
                versions: Dict[str, List[Tuple[Optional[packaging.version.Version], str]]] = {
                    name: [(None, version) for version in versions]
                    for name, versions in mod_dict.pop("versions").items()
                }
                intro = mod_dict.pop("intro")
                mod = BaseMultiqcModule(**mod_dict)
                mod.sections = sections
                mod.versions = versions
                mod.intro = intro
                log.info(f"Loading module {mod.name}")
                report.modules.append(mod)

            # Load data sources
            for mod_id, source_dict in data["report_data_sources"].items():
                for section, sources in source_dict.items():
                    for sname, source in sources.items():
                        report.data_sources[mod_id][section][sname] = source

            # Load plot_data dump
            for anchor, plot_dump in data["report_plot_data"].items():
                report.plot_data[anchor] = plot_dump

            # Load normalized plot inputs
            if "report_plot_input_data" in data and isinstance(data["report_plot_input_data"], dict):
                for anchor, rel_path in data["report_plot_input_data"].items():
                    prev_parquet_path = path.parent / rel_path
                    if not prev_parquet_path.exists():
                        log.error(f"Parquet file not found: {prev_parquet_path}")
                        continue
                    report.plot_input_data[anchor] = str(prev_parquet_path)

        except (json.JSONDecodeError, KeyError) as e:
            log.error(f"Error loading data from multiqc_data.json: {e}")
