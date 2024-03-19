import dataclasses
import io
import logging
import os
from typing import Dict, List, Union, Optional, Tuple

import math
import plotly.graph_objects as go

from multiqc.plots.plotly.plot import Plot, PlotType, BaseDataset
from multiqc.utils import util_functions, config
from multiqc.utils.config import update_dict

logger = logging.getLogger(__name__)


# {"name": "SAMPLE1", "color": "#111111", "data": [[x, y], [x, y], ...]}
LineT = Dict[str, Union[str, List[Tuple[Union[float, int, str], Union[float, int]]]]]


def plot(lists_of_lines: List[List[LineT]], pconfig: Dict) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param lists_of_lines: each dataset is a 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML with JS, ready to be inserted into the page
    """

    # if self.n_samples >= config.max_table_rows:
    # Get a line of median values at each point with an interval of max and
    # min values
    # Create a violin of median values in each sample, showing dots for outliers
    # Clicking on a dot of a violin will show the line plot for that sample

    p = LinePlot(pconfig, lists_of_lines)

    from multiqc.utils import report

    return p.add_to_report(report)


@dataclasses.dataclass
class Dataset(BaseDataset):
    lines: List[Dict]

    @staticmethod
    def create(
        dataset: BaseDataset,
        lines: List[Dict],
        pconfig: Dict,
    ) -> "Dataset":
        dataset: Dataset = Dataset(
            **dataset.__dict__,
            lines=lines,
        )

        # Prevent Plotly from parsing strings as numbers
        if pconfig.get("categories") or dataset.dconfig.get("categories"):
            dataset.layout["xaxis"]["type"] = "category"

        # convert HighCharts-style hardcoded trace parameters to Plotly style
        lines = []
        for src_line in dataset.lines:
            new_line = {
                "name": src_line["name"],
                "data": src_line["data"],
                "color": src_line.get("color"),
                "showlegend": src_line.get("showlegend", src_line.get("showInLegend", True)),
                "line": {
                    "dash": convert_dash_style(src_line.get("dash", src_line.get("dashStyle"))),
                    "width": src_line.get("line", {}).get("width", src_line.get("lineWidth")),
                },
            }
            if "marker" in src_line:
                new_line["mode"] = "lines+markers"
                new_line["marker"] = {
                    "symbol": src_line["marker"].get("symbol"),
                    "line": {
                        "width": src_line["marker"].get("line", {}).get("width", src_line["marker"].get("lineWidth")),
                        "color": src_line["marker"].get("line", {}).get("color", src_line["marker"].get("lineColor")),
                    },
                }
            lines.append(remove_nones_and_empty_dicts(new_line))

        dataset.lines = lines

        mode = pconfig.get("style", "lines")
        if config.lineplot_style == "lines+markers":
            mode = "lines+markers"

        dataset.trace_params.update(
            mode=mode,
            line={"width": 2},
        )
        if mode == "lines+markers":
            dataset.trace_params.update(
                line={"width": 0.6},
                marker={"size": 5},
            )
        return dataset

    def create_figure(
        self,
        layout: go.Layout,
        is_log=False,
        is_pct=False,
    ) -> go.Figure:
        """
        Create a Plotly figure for a dataset
        """
        layout.height += len(self.lines) * 5  # extra space for legend

        fig = go.Figure(layout=layout)
        for line in self.lines:
            xs = [x[0] for x in line["data"]]
            ys = [x[1] for x in line["data"]]
            params = dict(
                marker=line.get("marker", {}),
                line=line.get("line", {}),
                showlegend=line.get("showlegend", None),
                mode=line.get("mode", None),
            )
            params = update_dict(params, self.trace_params, none_only=True)
            params["marker"]["color"] = line.get("color")

            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    name=line["name"],
                    text=[line["name"]] * len(xs),
                    **params,
                )
            )
        return fig


class LinePlot(Plot):
    def __init__(self, pconfig: Dict, lists_of_lines: List[List[LineT]]):
        super().__init__(PlotType.LINE, pconfig, len(lists_of_lines))

        self.datasets: List[Dataset] = [
            Dataset.create(d, lines, pconfig) for d, lines in zip(self.datasets, lists_of_lines)
        ]

        # Make a tooltip always show on hover over any point on plot
        self.layout.hoverdistance = -1

        y_minrange = pconfig.get("y_minrange", pconfig.get("yMinRange"))
        x_minrange = pconfig.get("x_minrange", pconfig.get("xMinRange"))
        y_bands = pconfig.get("y_bands", pconfig.get("yPlotBands"))
        x_bands = pconfig.get("x_bands", pconfig.get("xPlotBands"))
        x_lines = pconfig.get("x_lines", pconfig.get("xPlotLines"))
        y_lines = pconfig.get("y_lines", pconfig.get("yPlotLines"))
        if y_minrange or y_bands or y_lines:
            # We don't want the bands to affect the calculated axis range, so we
            # find the min and the max from data points, and manually set the range.
            for dataset in self.datasets:
                minval = None
                maxval = None
                for line in dataset.lines:
                    ys = [x[1] for x in line["data"]]
                    if len(ys) > 0:
                        minval = min(ys) if minval is None else min(minval, min(ys))
                        maxval = max(ys) if maxval is None else max(maxval, max(ys))
                if maxval is not None and minval is not None:
                    maxval += (maxval - minval) * 0.05
                if minval is None:
                    minval = dataset.layout["yaxis"]["autorangeoptions"]["minallowed"]
                if maxval is None:
                    maxval = dataset.layout["yaxis"]["autorangeoptions"]["maxallowed"]
                clipmin = dataset.layout["yaxis"]["autorangeoptions"]["clipmin"]
                clipmax = dataset.layout["yaxis"]["autorangeoptions"]["clipmax"]
                if clipmin is not None and minval is not None and clipmin > minval:
                    minval = clipmin
                if clipmax is not None and maxval is not None and clipmax < maxval:
                    maxval = clipmax
                if y_minrange is not None and maxval is not None and minval is not None:
                    maxval = max(maxval, minval + y_minrange)
                if self.layout.yaxis.type == "log":
                    minval = math.log10(minval) if minval is not None and minval > 0 else None
                    maxval = math.log10(maxval) if maxval is not None and maxval > 0 else None
                dataset.layout["yaxis"]["range"] = [minval, maxval]

        if not pconfig.get("categories", False) and x_minrange or x_bands or x_lines:
            # same as above but for x-axis
            for dataset in self.datasets:
                minval = None
                maxval = None
                for line in dataset.lines:
                    xs = [x[0] for x in line["data"]]
                    if len(xs) > 0:
                        minval = min(xs) if minval is None else min(minval, min(xs))
                        maxval = max(xs) if maxval is None else max(maxval, max(xs))
                if minval is None:
                    minval = dataset.layout["xaxis"]["autorangeoptions"]["minallowed"]
                if maxval is None:
                    maxval = dataset.layout["xaxis"]["autorangeoptions"]["maxallowed"]
                clipmin = dataset.layout["xaxis"]["autorangeoptions"]["clipmin"]
                clipmax = dataset.layout["xaxis"]["autorangeoptions"]["clipmax"]
                if clipmin is not None and minval is not None and clipmin > minval:
                    minval = clipmin
                if clipmax is not None and maxval is not None and clipmax < maxval:
                    maxval = clipmax
                if x_minrange is not None and maxval is not None and minval is not None:
                    maxval = max(maxval, minval + x_minrange)
                if self.layout.xaxis.type == "log":
                    minval = math.log10(minval) if minval is not None and minval > 0 else None
                    maxval = math.log10(maxval) if maxval is not None and maxval > 0 else None
                dataset.layout["xaxis"]["range"] = [minval, maxval]

        self.layout.shapes = (
            [
                dict(
                    type="rect",
                    y0=band["from"],
                    y1=band["to"],
                    x0=0,
                    x1=1,
                    xref="paper",  # make x coords are relative to the plot paper [0,1]
                    fillcolor=band["color"],
                    line={
                        "width": 0,
                    },
                    layer="below",
                )
                for band in (y_bands or [])
            ]
            + [
                dict(
                    type="rect",
                    x0=band["from"],
                    x1=band["to"],
                    y0=0,
                    y1=1,
                    yref="paper",  # make y coords are relative to the plot paper [0,1]
                    fillcolor=band["color"],
                    line={
                        "width": 0,
                    },
                    layer="below",
                )
                for band in (x_bands or [])
            ]
            + [
                dict(
                    type="line",
                    xref="paper",
                    yref="y",
                    x0=0,
                    y0=line["value"],
                    x1=1,
                    y1=line["value"],
                    line={
                        "width": line.get("width", 2),
                        "dash": convert_dash_style(line.get("dash", line.get("dashStyle"))),
                        "color": line["color"],
                    },
                )
                for line in (y_lines or [])
            ]
            + [
                dict(
                    type="line",
                    yref="paper",
                    xref="x",
                    x0=line["value"],
                    y0=0,
                    x1=line["value"],
                    y1=1,
                    line={
                        "width": line.get("width", 2),
                        "dash": convert_dash_style(line.get("dash", line.get("dashStyle"))),
                        "color": line["color"],
                    },
                )
                for line in (x_lines or [])
            ]
        )

    @staticmethod
    def tt_label() -> Optional[str]:
        """Default tooltip label"""
        return "<br>%{x}: %{y}"

    def axis_controlled_by_switches(self) -> List[str]:
        """
        Return a list of axis names that are controlled by the log10 scale and percentage
        switch buttons
        """
        return ["yaxis"]

    def save_data_file(self, dataset: Dataset) -> None:
        y_by_x_by_sample = dict()
        last_cats = None
        shared_cats = True
        for line in dataset.lines:
            y_by_x_by_sample[line["name"]] = dict()

            # Check to see if all categories are the same
            if len(line["data"]) > 0 and isinstance(line["data"][0], list):
                if last_cats is None:
                    last_cats = [x[0] for x in line["data"]]
                elif last_cats != [x[0] for x in line["data"]]:
                    shared_cats = False

            for i, x in enumerate(line["data"]):
                if isinstance(x, list):
                    y_by_x_by_sample[line["name"]][x[0]] = x[1]
                else:
                    try:
                        y_by_x_by_sample[line["name"]][dataset.dconfig["categories"][i]] = x
                    except (ValueError, KeyError, IndexError):
                        y_by_x_by_sample[line["name"]][str(i)] = x

        # Custom tsv output if the x-axis varies
        if not shared_cats and config.data_format in ["tsv", "csv"]:
            sep = "\t" if config.data_format == "tsv" else ","
            fout = ""
            for line in dataset.lines:
                fout += line["name"] + sep + "X" + sep + sep.join([str(x[0]) for x in line["data"]]) + "\n"
                fout += line["name"] + sep + "Y" + sep + sep.join([str(x[1]) for x in line["data"]]) + "\n"

            fn = f"{dataset.uid}.{config.data_format_extensions[config.data_format]}"
            fpath = os.path.join(config.data_dir, fn)
            with io.open(fpath, "w", encoding="utf-8") as f:
                f.write(fout.encode("utf-8", "ignore").decode("utf-8"))
        else:
            util_functions.write_data_file(y_by_x_by_sample, dataset.uid)


def convert_dash_style(dash_style: str) -> str:
    """Convert dash style from Highcharts to Plotly"""
    mapping = {
        "Solid": "solid",
        "ShortDash": "dash",
        "ShortDot": "dot",
        "ShortDashDot": "dashdot",
        "ShortDashDotDot": "dashdot",
        "Dot": "dot",
        "Dash": "dash",
        "DashDot": "dashdot",
        "LongDash": "longdash",
        "LongDashDot": "longdashdot",
        "LongDashDotDot": "longdashdot",
    }
    if dash_style in mapping.values():  # Plotly style?
        return dash_style
    elif dash_style in mapping.keys():  # Highcharts style?
        return mapping[dash_style]
    return "solid"


def remove_nones_and_empty_dicts(d: Dict) -> Dict:
    """Remove None and empty dicts from a dict recursively."""
    if not isinstance(d, Dict):
        return d
    return {k: remove_nones_and_empty_dicts(v) for k, v in d.items() if v is not None and v != {}}
