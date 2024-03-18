from typing import Dict

import plotly.graph_objects as go
from multiqc.plots.plotly.plot import PlotType, BaseDataset, Plot
from multiqc.plots.plotly import line


def load_plot_from_json(dump: Dict, dataset_id: int) -> go.Figure:
    """
    Load a Plot object from a JSON dictionary
    """
    plot_type = PlotType(dump["plot_type"])
    dataset_dump = dump["datasets"][dataset_id]
    dataset = BaseDataset(
        plot_id=dump["id"],
        label=str(dataset_id + 1),
        uid=dump["id"],
        dconfig=dataset_dump["dconfig"],
        layout=dataset_dump["layout"],
        trace_params=dataset_dump["trace_params"],
        pct_range=dataset_dump["pct_range"],
    )
    layout = go.Layout(dump["layout"])  # make a copy
    layout.update(**dataset_dump["layout"])
    layout.width = layout.width or Plot.FLAT_PLOT_WIDTH

    if plot_type == PlotType.LINE:
        return line.Dataset(
            **dataset.__dict__,
            lines=dataset_dump["lines"],
        ).create_figure(
            layout=layout,
            is_log=dump["l_active"],
            is_pct=dump["p_active"],
        )

    else:
        raise ValueError(f"Plot type {plot_type} is unknown or unsupported")

    # if plot_type == PlotType.BAR:
    #     return bar.BarPlot.fig_from_dump(dump)
    # if plot_type == PlotType.BOX:
    #     return box.BoxPlot.fig_from_dump(dump)
    # if plot_type == PlotType.SCATTER:
    #     return scatter.ScatterPlot.fig_from_dump(dump)
    # if plot_type == PlotType.SCATTER:
    #     return heatmap.HeatmapPlot.fig_from_dump(dump)
    # if plot_type == PlotType.VIOLIN:
    #     return violin.ViolinPlot.fig_from_dump(dump)
    #
    # plot = Plot()
    # plot.id = dump["id"]
    # plot.layout = go.Layout(dump["layout"])
    # plot.pct_axis_update = dump["pct_axis_update"]
    # plot.axis_controlled_by_switches = dump["axis_controlled_by_switches"]
    # plot.l_active = dump["l_active"]
    # plot.p_active = dump["p_active"]
    # plot.pconfig = dump["config"]
    # plot.pconfig["square"] = dump["square"]
