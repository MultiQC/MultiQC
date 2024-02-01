import copy
import dataclasses
import io
import logging
import os
from typing import Dict, List, Union, Optional
import plotly.graph_objects as go


from multiqc.templates.plotly.plots.plot import Plot, PlotType, BaseDataset
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

    # if self.n_samples >= config.max_table_rows:
    # Get a line of median values at each point with an interval of max and
    # min values
    # Create a violin of median values in each sample, showing dots for outliers
    # Clicking on a dot of a violin will show the line plot for that sample

    p = LinePlot(pconfig, lists_of_lines)

    from multiqc.utils import report

    return p.add_to_report(report)


class LinePlot(Plot):
    @dataclasses.dataclass
    class Dataset(BaseDataset):
        lines: List[Dict]
        categories: List[str]

        @staticmethod
        def create(
            dataset: BaseDataset,
            lines: List[Dict],
            pconfig: Dict,
        ) -> "LinePlot.Dataset":
            dataset = LinePlot.Dataset(
                **dataset.__dict__,
                lines=lines,
                categories=pconfig.get("categories", []),
            )
            if dataset.categories:
                dataset.layout.xaxis.type = "category"
                # check that all lines have the same number of categories
                assert all(len(line["data"]) == len(dataset.categories) for line in dataset.lines), dataset.uid

            # Make a tooltip always show on hover over any point on plot
            dataset.layout.hoverdistance = -1

            dataset.trace_params.update(
                mode="lines" if config.lineplot_style == "lines" else "lines+markers",
                line={"width": 2 if config.lineplot_style == "lines" else 0.6},
                marker={"size": 4},
            )

            y_min_range = pconfig.get("yMinRange")
            y_bands = pconfig.get("yPlotBands")
            x_bands = pconfig.get("xPlotBands")
            x_lines = pconfig.get("xPlotLines")
            y_lines = pconfig.get("yPlotLines")
            if y_min_range or y_bands or x_bands or x_lines or y_lines:
                # We don't want the bands to affect the calculated y-axis range. So we
                # call `fig.full_figure_for_development` to force-calculate the figure size
                # before we add the bands, and then force set the calculated range.
                dev_fig = dataset.create_figure(dataset.layout).full_figure_for_development(warn=False)
                dataset.layout.yaxis.range = dev_fig.layout.yaxis.range
                dataset.layout.xaxis.range = dev_fig.layout.xaxis.range
                if y_min_range:
                    if dataset.layout.yaxis.range[1] - dataset.layout.yaxis.range[0] < y_min_range:
                        dataset.layout.yaxis.range = [
                            dataset.layout.yaxis.range[0],
                            dataset.layout.yaxis.range[0] + y_min_range,
                        ]

                dataset.layout.shapes = (
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
                                "width": line["width"],
                                "dash": convert_dash_style(line.get("dashStyle")),
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
                                "width": line["width"],
                                "dash": convert_dash_style(line.get("dashStyle")),
                                "color": line["color"],
                            },
                        )
                        for line in (x_lines or [])
                    ]
                )

            return dataset

        def create_figure(
            self,
            layout: Optional[go.Layout] = None,
            is_log=False,
            is_pct=False,
        ) -> go.Figure:
            """
            Create a Plotly figure for a dataset
            """
            fig = go.Figure(layout=layout or self.layout)

            for line in self.lines:
                if len(line["data"]) > 0 and isinstance(line["data"][0], list):
                    xs = [x[0] for x in line["data"]]
                    ys = [x[1] for x in line["data"]]
                elif self.categories:
                    assert len(line["data"]) == len(self.categories)
                    xs = self.categories
                    ys = line["data"]
                else:
                    xs = [x for x in range(len(line["data"]))]
                    ys = line["data"]

                params = copy.deepcopy(self.trace_params)
                if "dashStyle" in line:
                    params["line"]["dash"] = convert_dash_style(line["dashStyle"].lower())
                if "lineWidth" in line:
                    params["line"]["width"] = line["lineWidth"]
                fig.add_trace(
                    go.Scatter(
                        x=xs,
                        y=ys,
                        name=line["name"],
                        text=[line["name"]] * len(xs),
                        marker_color=line.get("color"),
                        **params,
                    )
                )
            return fig

    def __init__(self, pconfig: Dict, lists_of_lines: List[List[LineT]]):
        super().__init__(PlotType.LINE, pconfig, len(lists_of_lines))

        self.datasets: List[LinePlot.Dataset] = [
            LinePlot.Dataset.create(d, lines, pconfig) for d, lines in zip(self.datasets, lists_of_lines)
        ]

        for dataset in self.datasets:
            if self.flat:
                dataset.layout.height += len(dataset.lines) * 5  # extra space for legend

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
                        y_by_x_by_sample[line["name"]][dataset.categories[i]] = x
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
    return {
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
    }.get(dash_style, "solid")
