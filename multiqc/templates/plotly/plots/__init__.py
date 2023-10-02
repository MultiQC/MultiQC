from typing import Dict, cast

import plotly.graph_objects as go

from multiqc.utils import mqc_colour


def basic_figure(pconfig: Dict, plt_height: int) -> go.Figure:
    """
    Plotly Figure object suitable for all plots.
    """
    layout = go.Layout(
        title=dict(
            text=pconfig["title"],
            xanchor="center",
            x=0.5,
            font=dict(size=20),
        ),
        yaxis=dict(
            title=dict(text=pconfig.get("ylab")),
            showgrid=False,
            categoryorder="category descending",
        ),
        xaxis=dict(
            title=dict(text=pconfig.get("xlab")),
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
        height=plt_height,
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
