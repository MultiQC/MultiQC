from dataclasses import dataclass
from typing import Dict, cast

import plotly.graph_objects as go

from multiqc.utils import config, mqc_colour, report

# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod = None


def get_template_mod():
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    return _template_mod


@dataclass
class PlotSettings:
    def __init__(self, pconfig: Dict):
        self.title = pconfig.get("title")
        self._id = pconfig.get("id")
        self.xlab = pconfig.get("xlab")
        self.ylab = pconfig.get("ylab")
        self.height = pconfig.get("height")

        self.add_pct_tab = pconfig.get("cpswitch", True) is not False
        self.add_log_tab = pconfig.get("logswitch", False)

        # Counts / Percentages Switch
        self.c_label = pconfig.get("cpswitch_counts_label", "Counts")
        self.p_label = pconfig.get("cpswitch_percent_label", "Percentages")
        self.l_label = pconfig.get("logswitch_label", "Log10")
        self.c_active = ""
        self.p_active = ""
        self.l_active = ""
        self.stacking = pconfig.get("stacking", "normal")
        if self.add_pct_tab and not config.simple_output:
            if pconfig.get("logswitch_active") is True:
                self.l_active = "active"
            elif pconfig.get("cpswitch_c_active", True) is True:
                self.c_active = "active"
            else:
                self.p_active = "active"
                self.stacking = "percent"

        self.data_labels = pconfig.get("data_labels", [])
        self.save_data_file = pconfig.get("save_data_file", True)

    @property
    def id(self) -> str:
        return self._id

    @id.setter
    def id(self, id: str):
        # Sanitise plot ID and check for duplicates
        self._id = report.save_htmlid(id)


def basic_figure(settings: PlotSettings) -> go.Figure:
    """
    Plotly Figure object suitable for all plots.
    """
    layout = go.Layout(
        title=dict(
            text=settings.title,
            xanchor="center",
            x=0.5,
            font=dict(size=20),
        ),
        yaxis=dict(
            title=dict(text=settings.ylab),
            showgrid=False,
            categoryorder="category descending",
        ),
        xaxis=dict(
            title=dict(text=settings.xlab),
            gridcolor="rgba(0,0,0,0.1)",
        ),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        hovermode="y unified",
        # template="simple_white",
        font={"color": "Black", "family": "Lucida Grande"},
        colorway=mqc_colour.mqc_colour_scale.COLORBREWER_SCALES["plot_defaults"],
    )
    fig = go.Figure()
    fig.update_layout(
        # The function expects a dict, even though go.Layout works just fine
        cast(dict, layout),
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.2,
            xanchor="center",
            x=0.5,
        ),
        barmode="stack",
        height=settings.height,
    )

    fig.add_annotation(
        xanchor="right",
        yanchor="bottom",
        x=1.1,
        y=-0.3,
        text="Created with MultiQC",
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(size=10, color="rgba(0,0,0,0.5)"),
    )

    return fig
