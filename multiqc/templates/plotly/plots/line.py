import base64
import io
import logging
import os
from pathlib import Path
from typing import Dict, List, Union, Optional

import math
import plotly.graph_objects as go

from multiqc.templates.plotly.plots import get_template_mod
from multiqc.templates.plotly.plots.plot import View, Plot
from multiqc.utils import util_functions, config

logger = logging.getLogger(__name__)


# {"name": "SAMPLE1", "color": "#111111", "data": [[x, y], [x, y], ...]}
SampleLineT = Dict[str, Union[str, List]]


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

    def _save_data_file(
        self,
        pid: str,
        dataset: List[SampleLineT],
    ) -> None:
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
            with io.open(os.path.join(config.data_dir, f"{pid}.txt"), "w", encoding="utf-8") as f:
                print(fout.encode("utf-8", "ignore").decode("utf-8"), file=f)
        else:
            util_functions.write_data_file(fdata, pid)

    def _flat_imgs_for_dataset(
        self,
        pidx: int,
        pid: str,
        dataset: List[SampleLineT],
    ) -> str:
        """
        Build a static images for different views of a dataset (counts, log scale),
        return an HTML wrapper.
        """
        # Save plot data to file
        if self.save_data_file:
            self._save_data_file(pid, dataset)

        # Calculate log10 values
        if self.add_log_tab:
            for ds in dataset:
                ds["data_log"] = []
                for x, y in ds["data"]:
                    if y == 0:
                        y = 0
                    else:
                        y = math.log10(y)
                    ds["data_log"].append([x, y])

        views = [
            View(
                dataset,
                active=not self.l_active,
                suffix="",
                label=self.c_label,
                xaxis_tickformat="",
            ),
        ]
        if self.add_log_tab:
            views.append(
                View(
                    [
                        {
                            "data": ds["data_log"],
                            "name": ds["name"],
                            "color": ds["color"],
                        }
                        for ds in dataset
                    ],
                    active=False,
                    suffix="_log",
                    label=self.l_label,
                    xaxis_tickformat="",
                )
            )

        html = ""
        for view in views:
            html += self._flat_img_for_view(
                view,
                pidx,
                f"{pid}{view.suffix}",
            )
        return html

    def _flat_img_for_view(
        self,
        view: View,
        pidx: int,
        pid: str,
    ) -> str:
        """
        Build one static image, return an HTML wrapper.
        """
        pid = f"{pid}{view.suffix}"

        # Should this plot be hidden on report load?
        hide_div = ""
        if pidx > 0 or not view.active:
            hide_div = ' style="display:none;"'

        fig = go.Figure(layout=self.layout())
        for sdata in view.data:
            if len(sdata["data"]) > 0 and isinstance(sdata["data"][0], list):
                x = [x[0] for x in sdata["data"]]
                y = [x[1] for x in sdata["data"]]
            else:
                x = [x for x in range(len(sdata["data"]))]
                y = sdata["data"]

            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    name=sdata["name"],
                    mode="lines+markers",
                    marker=dict(size=5),
                )
            )

        # Save the plot to the data directory if export is requested
        if config.export_plots:
            for fformat in config.export_plot_formats:
                # Make the directory if it doesn't already exist
                plot_dir = os.path.join(config.plots_dir, fformat)
                if not os.path.exists(plot_dir):
                    os.makedirs(plot_dir)
                # Save the plot
                plot_fn = os.path.join(plot_dir, f"{pid}.{fformat}")
                fig.write_image(
                    plot_fn,
                    format=fformat,
                    width=fig.layout.width,
                    height=fig.layout.height,
                    scale=1,
                )

        # Output the figure to a base64 encoded string
        if getattr(get_template_mod(), "base64_plots", True) is True:
            img_buffer = io.BytesIO()
            fig.write_image(
                img_buffer,
                format="png",
                width=fig.layout.width,
                height=fig.layout.height,
                scale=1,
            )
            b64_img = base64.b64encode(img_buffer.getvalue()).decode("utf8")
            img_buffer.close()
            return f'<div class="mqc_mplplot" id="{pid}"{hide_div}><img src="data:image/png;base64,{b64_img}" /></div>'

        # Link to the saved image
        else:
            plot_relpath = Path(config.plots_dir_name) / "png" / f"{pid}.png"
            plot_relpath.parent.mkdir(parents=True, exist_ok=True)
            fig.write_image(
                plot_relpath,
                format="png",
                width=900,
                height=fig.layout.height,
                scale=1,
            )
            return f'<div class="mqc_mplplot" id="{pid}"{hide_div}><img src="{plot_relpath}" /></div>'
