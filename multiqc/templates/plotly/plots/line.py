import dataclasses
import io
import logging
import os
from typing import Dict, List, Union

import math
import plotly.graph_objects as go

from multiqc.templates.plotly.plots.plot import Plot, PlotType, Dataset
from multiqc.utils import util_functions, config

logger = logging.getLogger(__name__)


# {"name": "SAMPLE1", "color": "#111111", "data": [[x, y], [x, y], ...]}
LineT = Dict[str, Union[str, List[List[float]]]]


def plot(lists_of_lines: List[List[LineT]], pconfig: Dict) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param lists_of_lines: each dataset is a 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML with JS, ready to be inserted into the page
    """
    p = LinePlot(pconfig, lists_of_lines)

    from multiqc.utils import report

    return p.add_to_report(report)


@dataclasses.dataclass
class LineDataset(Dataset):
    """Bar dataset should also carry the list of samples"""

    lines: List[Dict]


class LinePlot(Plot):
    def __init__(self, pconfig: Dict, lists_of_lines: List[List[LineT]]):
        super().__init__(PlotType.LINE, pconfig, len(lists_of_lines))

        # Extend each dataset object with a list of samples
        self.datasets: List[LineDataset] = [
            LineDataset(**d.__dict__, lines=lines) for d, lines in zip(self.datasets, lists_of_lines)
        ]

        self.categories: List[str] = pconfig.get("categories", [])

        # Make a tooltip always show on hover over any point on plot
        self.layout.hoverdistance = -1
        # A tooltip will show numbers for all lines crossing this vertical line
        # self.layout.hovermode = "x"
        # Default precision for floating numbers is too high - allowing to override it
        if "tt_decimals" in pconfig:
            self.layout.yaxis.hoverformat = f".{pconfig['tt_decimals']}f"
        # self.tt_suffix: str = pconfig.get("tt_suffix", "")
        # self.tt_label: str = pconfig.get(
        #     "tt_label",
        #     f"%{{x}}: %{{y:,.{self.tt_decimals}f}}{self.tt_suffix}"
        # )
        if self.flat and self.datasets:
            self.layout.height += len(self.datasets[0].lines) * 5  # extra space for legend

        y_bands = self.pconfig.get("yPlotBands")
        x_bands = self.pconfig.get("xPlotBands")
        if y_bands or x_bands:
            # We don't want the bands to affect the calculated y-axis range. So we
            # call `fig.full_figure_for_development` to force-calculate the figure size
            # before we add the bands, and then force set the calculated range.
            dev_figs = [self.create_figure(self.layout, d) for d in self.datasets]
            dev_figs = [fig.full_figure_for_development(warn=False) for fig in dev_figs]
            self.layout.yaxis.range = [
                min([fig.layout.yaxis.range[0] for fig in dev_figs]),
                max([fig.layout.yaxis.range[1] for fig in dev_figs]),
            ]
            self.layout.xaxis.range = [
                min([fig.layout.xaxis.range[0] for fig in dev_figs]),
                max([fig.layout.xaxis.range[1] for fig in dev_figs]),
            ]

            self.layout.shapes = [
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
                for band in self.pconfig.get("yPlotBands", [])
            ] + [
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
                for band in self.pconfig.get("xPlotBands", [])
            ]

        self.trace_params = dict(
            mode="lines+markers",
            line=dict(width=0.6),
            marker=dict(size=4),
        )
        if config.lineplot_style == "lines":
            self.trace_params = dict(
                mode="lines",
            )

    def dump_for_javascript(self):
        """Serialise the data to pick up in plotly-js"""
        d = super().dump_for_javascript()
        d["categories"] = self.categories
        return d

    def create_figure(self, layout: go.Layout, dataset: LineDataset, is_log=False, is_pct=False):
        """
        Create a Plotly figure for a dataset
        """

        import json

        with open(f"/Users/vlad/git/playground/dumps/{self.id}-layout.json", "w") as f:
            f.write(json.dumps(layout.to_plotly_json()))
        with open(f"/Users/vlad/git/playground/dumps/{self.id}-lines.json", "w") as f:
            f.write(json.dumps(dataset.lines))
        with open(f"/Users/vlad/git/playground/dumps/{self.id}-pconfig.json", "w") as f:
            f.write(json.dumps(self.pconfig))

        if self.id == "fastqc_per_base_n_content_plot":
            print()

        fig = go.Figure(layout=layout)
        for line in dataset.lines:
            if len(line["data"]) > 0 and isinstance(line["data"][0], list):
                xs = [x[0] for x in line["data"]]
                ys = [x[1] for x in line["data"]]
            else:
                xs = [x for x in range(len(line["data"]))]
                ys = line["data"]
            if is_log:
                ys = [math.log10(y) if y else 0 for y in ys]

            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    name=line["name"],
                    marker_color=line.get("color"),
                    **self.trace_params,
                    # mode="lines+markers",
                    # line=dict(
                    #     width=0.6,
                    # ),
                    # marker=dict(
                    #     size=4,
                    #     color=line.get("color"),
                    # ),
                )
            )
        return fig

    def save_data_file(self, dataset: LineDataset) -> None:
        fdata = dict()
        last_cats = None
        shared_cats = True
        for line in dataset.lines:
            fdata[line["name"]] = dict()

            # Check to see if all categories are the same
            if len(line["data"]) > 0 and isinstance(line["data"][0], list):
                if last_cats is None:
                    last_cats = [x[0] for x in line["data"]]
                elif last_cats != [x[0] for x in line["data"]]:
                    shared_cats = False

            for i, x in enumerate(line["data"]):
                if isinstance(x, list):
                    fdata[line["name"]][x[0]] = x[1]
                else:
                    try:
                        fdata[line["name"]][self.categories[i]] = x
                    except (ValueError, KeyError, IndexError):
                        fdata[line["name"]][str(i)] = x

        # Custom tsv output if the x-axis varies
        if not shared_cats and config.data_format == "tsv":
            fout = ""
            for line in dataset.lines:
                fout += "\t" + "\t".join([str(x[0]) for x in line["data"]])
                fout += "\n{}\t".format(line["name"])
                fout += "\t".join([str(x[1]) for x in line["data"]])
                fout += "\n"
            with io.open(os.path.join(config.data_dir, f"{dataset.uid}.txt"), "w", encoding="utf-8") as f:
                f.write(fout.encode("utf-8", "ignore").decode("utf-8"))
        else:
            util_functions.write_data_file(fdata, dataset.uid)
