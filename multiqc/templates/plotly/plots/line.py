import io
import logging
import os
from typing import Dict, List, Union, Optional

import math
import plotly.graph_objects as go

from multiqc.templates.plotly.plots.plot import Plot
from multiqc.utils import util_functions, config

logger = logging.getLogger(__name__)


# {"name": "SAMPLE1", "color": "#111111", "data": [[x, y], [x, y], ...]}
SampleLineT = Dict[str, Union[str, List[List[float]]]]


def plot(
    datasets: List[List[SampleLineT]],
    pconfig: Dict,
) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param datasets: each dataset is a 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML with JS, ready to be inserted into the page
    """
    from multiqc.utils import report

    p = LinePlot(pconfig, len(datasets))

    return p.add_to_report(datasets, report)


class LinePlot(Plot):
    def __init__(self, pconfig: Dict, *args):
        super().__init__("xy_line", pconfig, *args)

        self.categories: List[str] = pconfig.get("categories", [])

        self.tt_decimals: Optional[int] = pconfig.get("tt_decimals")
        self.tt_suffix: str = pconfig.get("tt_suffix", "")
        # self.tt_label: str = pconfig.get(
        #     "tt_label",
        #     f"%{{x}}: %{{y:,.{self.tt_decimals}f}}{self.tt_suffix}"
        # )

    def layout(self) -> go.Layout:
        layout: go.Layout = super().layout()

        layout.showlegend = False

        # Make a tooltip always show on hover over any point on plot
        layout.hoverdistance = -1
        # A tooltip will show numbers for all lines crossing this vertical line
        layout.hovermode = "x"
        # Default precision for floating numbers is too high - allowing to override it
        if self.tt_decimals is not None:
            layout.yaxis.hoverformat = f".{self.tt_decimals}f"

        return layout

    def populate_figure(self, fig: go.Figure, dataset: List[SampleLineT], is_log=False, is_pct=False):
        for sample in dataset:
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

    def save_data_file(self, dataset: List[SampleLineT], uid: str) -> None:
        fdata = dict()
        last_cats = None
        shared_cats = True
        for ds in dataset:
            fdata[ds["name"]] = dict()

            # Check to see if all categories are the same
            if len(ds["data"]) > 0 and isinstance(ds["data"][0], list):
                if last_cats is None:
                    last_cats = [x[0] for x in ds["data"]]
                elif last_cats != [x[0] for x in ds["data"]]:
                    shared_cats = False

            for i, x in enumerate(ds["data"]):
                if isinstance(x, list):
                    fdata[ds["name"]][x[0]] = x[1]
                else:
                    try:
                        fdata[ds["name"]][self.categories[i]] = x
                    except Exception:
                        fdata[ds["name"]][str(i)] = x

        # Custom tsv output if the x-axis varies
        if not shared_cats and config.data_format == "tsv":
            fout = ""
            for ds in dataset:
                fout += "\t" + "\t".join([str(x[0]) for x in ds["data"]])
                fout += "\n{}\t".format(ds["name"])
                fout += "\t".join([str(x[1]) for x in ds["data"]])
                fout += "\n"
            with io.open(os.path.join(config.data_dir, f"{uid}.txt"), "w", encoding="utf-8") as f:
                f.write(fout.encode("utf-8", "ignore").decode("utf-8"))
        else:
            util_functions.write_data_file(fdata, uid)
