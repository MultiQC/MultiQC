"""Plotly bargraph functionality."""
# import base64

import plotly.graph_objects as go


def plotly_bargraph(plot_data, plot_series, pconfig):
    # print(plot_data)
    # print(plot_series)
    # print(pconfig)

    # NOTE: Only first plot for now, if multiple plots
    # Height has a default, then adjusted by the number of samples
    plt_height = len(plot_series[0]) / 186  # Default, empirically determined
    plt_height = max(600, plt_height)  # At least 512px tall
    plt_height = min(2560, plt_height)  # Cap at 2560 tall
    fig = go.Figure()
    for data in plot_data[0]:
        fig.add_trace(
            go.Bar(
                y=plot_series[0],
                x=data["data"],
                name=data["name"],
                orientation="h",
                marker=dict(
                    color=data.get("color"),
                    line=dict(width=0),
                ),
            )
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
    fig.update_layout(
        go.Layout(
            title=dict(
                text=pconfig["title"],
                xanchor="center",
                x=0.5,
                font=dict(size=20),
            ),
            yaxis=dict(
                title=dict(text=pconfig.get("xlab")),
                showgrid=False,
                categoryorder="category descending",
            ),
            xaxis=dict(
                title=dict(text=pconfig.get("ylab")),
                gridcolor="rgba(0,0,0,0.1)",
            ),
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
            hovermode="y unified",
            # template="simple_white",
            font={"color": "Black", "family": "Lucida Grande"},
            colorway=[
                "#7cb5ec",
                "#434348",
                "#90ed7d",
                "#f7a35c",
                "#8085e9",
                "#f15c80",
                "#e4d354",
                "#2b908f",
                "#f45b5b",
                "#91e8e1",
            ],
        ),
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.2,
            xanchor="center",
            x=0.5,
        ),
        barmode="stack",
        # width=1200,
        height=plt_height,
    )

    # img_bytes = fig.to_image(format="png", scale=2)
    # img_uri = "data:image/png;base64," + base64.b64encode(img_bytes).decode("utf8")
    # return f'<div class="mqc_mplplot"><img src="{img_uri}" /></div>'

    img_html = fig.to_html(full_html=False, include_plotlyjs="cdn")
    return f'<div class="mqc_mplplot">{img_html}</div>'
