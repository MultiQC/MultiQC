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
SampleLineT = Dict[str, Union[str, List[List[float]]]]


def plot(datasets: List[List[SampleLineT]], pconfig: Dict) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param datasets: each dataset is a 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML with JS, ready to be inserted into the page
    """
    p = LinePlot(pconfig, datasets)

    from multiqc.utils import report

    return p.add_to_report(report)


class LinePlot(Plot):
    def __init__(self, pconfig: Dict, datasets: List):
        super().__init__(PlotType.LINE, pconfig, datasets)

        self.categories: List[str] = pconfig.get("categories", [])

        # Make a tooltip always show on hover over any point on plot
        self.layout.hoverdistance = -1
        # A tooltip will show numbers for all lines crossing this vertical line
        self.layout.hovermode = "x"
        # Default precision for floating numbers is too high - allowing to override it
        if "tt_decimals" in pconfig:
            self.layout.yaxis.hoverformat = f".{pconfig['tt_decimals']}f"
        # self.tt_suffix: str = pconfig.get("tt_suffix", "")
        # self.tt_label: str = pconfig.get(
        #     "tt_label",
        #     f"%{{x}}: %{{y:,.{self.tt_decimals}f}}{self.tt_suffix}"
        # )
        if self.flat and self.datasets:
            self.layout.height += len(self.datasets[0].data) * 5  # extra space for legend

    def create_figure(self, layout: go.Layout, dataset: Dataset, is_log=False, is_pct=False):
        """
        Create a Plotly figure for a dataset
        """

        # import json
        #
        # with open(f"/Users/vlad/git/playground/{self.id}-layout.json", "w") as f:
        #     f.write(json.dumps(layout.to_plotly_json()))
        # with open(f"/Users/vlad/git/playground/{self.id}-data.json", "w") as f:
        #     f.write(json.dumps(dataset.data))

        fig = go.Figure(layout=layout)
        data: List[SampleLineT] = dataset.data
        for sample in data:
            if len(sample["data"]) > 0 and isinstance(sample["data"][0], list):
                xs = [x[0] for x in sample["data"]]
                ys = [x[1] for x in sample["data"]]
            else:
                xs = [x for x in range(len(sample["data"]))]
                ys = sample["data"]
            if is_log:
                ys = [math.log10(y) if y else 0 for y in ys]

            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    name=sample["name"],
                    mode="lines+markers",
                    marker=dict(size=5),
                )
            )
        return fig

    def save_data_file(self, dataset: Dataset) -> None:
        fdata = dict()
        last_cats = None
        shared_cats = True
        for sample in dataset.data:
            fdata[sample["name"]] = dict()

            # Check to see if all categories are the same
            if len(sample["data"]) > 0 and isinstance(sample["data"][0], list):
                if last_cats is None:
                    last_cats = [x[0] for x in sample["data"]]
                elif last_cats != [x[0] for x in sample["data"]]:
                    shared_cats = False

            for i, x in enumerate(sample["data"]):
                if isinstance(x, list):
                    fdata[sample["name"]][x[0]] = x[1]
                else:
                    try:
                        fdata[sample["name"]][self.categories[i]] = x
                    except Exception:
                        fdata[sample["name"]][str(i)] = x

        # Custom tsv output if the x-axis varies
        if not shared_cats and config.data_format == "tsv":
            fout = ""
            for sample in dataset.data:
                fout += "\t" + "\t".join([str(x[0]) for x in sample["data"]])
                fout += "\n{}\t".format(sample["name"])
                fout += "\t".join([str(x[1]) for x in sample["data"]])
                fout += "\n"
            with io.open(os.path.join(config.data_dir, f"{dataset.uid}.txt"), "w", encoding="utf-8") as f:
                f.write(fout.encode("utf-8", "ignore").decode("utf-8"))
        else:
            util_functions.write_data_file(fdata, dataset.uid)
