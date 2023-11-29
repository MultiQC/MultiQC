"""Plotly bargraph functionality."""
import base64
import io
import logging
import os
from pathlib import Path
from typing import Dict, List

import math
import plotly.graph_objects as go

from multiqc.templates.plotly.plots import get_template_mod
from multiqc.templates.plotly.plots.plot import View, Plot
from multiqc.utils import config, util_functions

logger = logging.getLogger(__name__)


def plot(
    datasets: List[List[Dict]],
    samples_lists: List[List[str]],
    pconfig: Dict,
) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param datasets: each dataset is a list of dicts with the keys: {name, color, data},
        where `name` is the category name, `color` is the color of the bar,
        and `data` is a list of values for each sample. Each outer list will
        correspond a separate tab.
    :param samples_lists: list of lists of bar names (that is, sample names). Similarly,
        each outer list will correspond to a separate tab.
    :param pconfig: Plot configuration dictionary
    :return: HTML with JS, ready to be inserted into the page
    """
    from multiqc.utils import report

    for samples, ds in zip(samples_lists, datasets):
        for cat in ds:
            cat["samples"] = samples

    p = BarPlot(pconfig, len(datasets), max(len(samples) for samples in samples_lists))

    # Calculate and save percentages
    if p.add_pct_tab:
        for pidx, data_by_cat in enumerate(datasets):
            # Count totals for each category
            sums = [0 for _ in data_by_cat[0]["data"]]
            for cat_idx, d in enumerate(data_by_cat):
                for sample_idx, v in enumerate(d["data"]):
                    if not math.isnan(v):
                        sums[sample_idx] += v
            # Now, calculate percentages for each category
            for cat_idx, d in enumerate(data_by_cat):
                values = [x for x in d["data"]]
                if len(values) < len(d["samples"]):
                    values.extend([0] * (len(d["samples"]) - len(values)))
                for key, var in enumerate(values):
                    sum_for_cat = sums[key]
                    if sum_for_cat == 0:
                        values[key] = 0
                    else:
                        values[key] = float(var + 0.0) / float(sum_for_cat)
                d["data_pct"] = values

    return p.add_to_report(datasets, report)


class BarPlot(Plot):
    def __init__(self, pconfig: Dict, max_n_samples: int, *args):
        super().__init__("bar_graph", pconfig, *args)

        self.add_pct_tab = pconfig.get("cpswitch", True) is not False
        self.p_label = pconfig.get("cpswitch_percent_label", "Percentages")
        self.p_active = ""
        self.stacking = pconfig.get("stacking", "stack")
        if self.add_pct_tab and pconfig.get("cpswitch_c_active", True) is not True and not config.simple_output:
            self.c_active = ""
            self.p_active = "active"
            self.stacking = "relative"

        if not self.height:
            # Height has a default, then adjusted by the number of samples
            self.height = max_n_samples // 186  # Default, empirically determined
            self.height = max(600, self.height)
            self.height = min(2560, self.height)

    def layout(self) -> go.Layout:
        layout: go.Layout = super().layout()
        layout.update(
            {
                "barmode": self.stacking,
                "hovermode": "y unified",
                "showlegend": True,
                "yaxis": dict(
                    # the plot is "transposed", so yaxis corresponds to the horizontal axis
                    title=dict(text=self.xlab),
                    showgrid=False,
                    categoryorder="category descending",
                    automargin=True,
                ),
                "xaxis": dict(
                    title=dict(text=self.ylab),
                ),
                "legend": dict(
                    orientation="h",
                    yanchor="top",
                    y=-0.2,
                    xanchor="center",
                    x=0.5,
                ),
            }
        )
        return layout

    # def flat_plot(
    #     self,
    #     datasets: List,
    #     pid: str,
    #     pids: List[str],
    # ) -> str:
    #     """
    #     Create an HTML for a static plot version.
    #     """
    #     html = super().flat_plot(datasets, pid, pids)
    #
    #     # Buttons to cycle through different datasets
    #     if len(datasets) > 1 and not config.simple_output:
    #         # Collect all categories, and fill in with zeroes for samples that having
    #         # any of cats missing
    #         cat_to_color = dict()
    #         for ds in datasets:
    #             for d in ds:
    #                 cat_to_color[d["name"]] = d["color"]
    #         for data_by_cat in datasets:
    #             for cat, color in cat_to_color.items():
    #                 samples = data_by_cat[0]["samples"]
    #                 if cat not in [d["name"] for d in data_by_cat]:
    #                     data_by_cat.append(
    #                         {
    #                             "name": cat,
    #                             "color": color,
    #                             "data": [0 for _ in samples],
    #                             "data_pct": [0.0 for _ in samples],
    #                             "data_log": [0.0 for _ in samples],
    #                         }
    #                     )
    #         # Sort categories by name
    #         for data_by_cat in datasets:
    #             data_by_cat.sort(key=lambda x: x["name"])
    #
    #     # Finally, build and save plots
    #     for pidx, (pid, data_by_cat) in enumerate(zip(pids, datasets)):
    #         html += self.dataset_to_imgs(pidx, pid, data_by_cat)
    #
    #     # Close wrapping div
    #     html += "</div>"
    #     return html

    @staticmethod
    def _save_data_file(data_by_cat, pid):
        fdata = {}
        for d in data_by_cat:
            for d_idx, d_val in enumerate(d["data"]):
                s_name = d["samples"][d_idx]
                if s_name not in fdata:
                    fdata[s_name] = dict()
                fdata[s_name][d["name"]] = d_val
        util_functions.write_data_file(fdata, pid)

    def _flat_imgs_for_dataset(self, pidx, pid, data_by_cat) -> str:
        """
        Build a static images for different views of a dataset (counts, percentages, log scale),
        return an HTML wrapper.
        """
        # Save plot data to file
        if self.save_data_file:
            self._save_data_file(data_by_cat, pid)

        # # Switch out NaN for 0s
        # for idx, d in enumerate(data_by_cat):
        #     data_by_cat[idx]["data"] = [x if not math.isnan(x) else 0 for x in d["data"]]

        # Calculate log10 values
        if self.add_log_tab:
            for cat_idx, d in enumerate(data_by_cat):
                values = [x for x in d["data"]]
                if len(values) < len(d["samples"]):
                    values.extend([0] * (len(d["samples"]) - len(values)))
                for key, var in enumerate(values):
                    if var == 0:
                        values[key] = 0
                    else:
                        values[key] = math.log10(var)
                d["data_log"] = values

        views = [
            View(
                data_by_cat,
                active=not self.p_active,
                suffix="",
                label=self.c_label,
                xaxis_tickformat="",
            ),
        ]
        if self.add_pct_tab:
            views.append(
                View(
                    [
                        {
                            "data": d["data_pct"],
                            "samples": d["samples"],
                            "name": d["name"],
                            "color": d["color"],
                        }
                        for d in data_by_cat
                    ],
                    active=self.p_active,
                    suffix="_pc",
                    label=self.p_label,
                    xaxis_tickformat=".0%",
                )
            )
        if self.add_log_tab:
            views.append(
                View(
                    [
                        {
                            "data": d["data_log"],
                            "samples": d["samples"],
                            "name": d["name"],
                            "color": d["color"],
                        }
                        for d in data_by_cat
                    ],
                    active=False,
                    suffix="_log",
                    label=self.l_label,
                    xaxis_tickformat="",
                )
            )

        html = ""
        for view in views:
            html += self._view_to_img(
                view,
                pidx,
                f"{pid}{view.suffix}",
            )
        return html

    def _view_to_img(
        self,
        view: View,
        pidx: int,
        pid: str,
    ) -> str:
        """
        Build one static image, return an HTML wrapper.
        """
        # Should this plot be hidden on report load?
        hide_div = ""
        if pidx > 0 or not view.active:
            hide_div = ' style="display:none;"'

        fig = go.Figure(layout=self.layout())
        for d in view.data:
            fig.add_trace(
                go.Bar(
                    y=d["samples"],
                    x=d["data"],
                    name=d["name"],
                    orientation="h",
                    marker=dict(color=d.get("color"), line=dict(width=0)),
                ),
            )
        if view.xaxis_tickformat:
            fig.update_layout({"xaxis": {"tickformat": view.xaxis_tickformat}})

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
