"""
ECharts plotting backend: Python option builders.

All ECharts knowledge lives in this one package, one module per plot type,
dispatched on `PlotType`. See "The ECharts model->JSON contract" in
`multiqc-echarts-exploration/BUILD_PLAN.md` for the full serialize/get_option
contract this package implements.

Circular-import note: `multiqc/plots/plot.py` must import this package LAZILY
(inside function bodies only), because `bar.py` imports `multiqc.plots.bargraph`,
which imports `multiqc.plots.plot`.
"""

import json
from types import ModuleType
from typing import TYPE_CHECKING, Any, Dict

from multiqc import config
from multiqc.plots.echarts import bar, box, heatmap, line, scatter
from multiqc.types import PlotType

if TYPE_CHECKING:
    from multiqc.plots.plot import Plot

_BUILDERS: Dict[PlotType, ModuleType] = {
    PlotType.BAR: bar,
    PlotType.LINE: line,
    PlotType.SCATTER: scatter,
    PlotType.HEATMAP: heatmap,
    PlotType.BOX: box,
}


def _builder(plot_type: PlotType) -> ModuleType:
    # `Plot.plot_type` is declared as `PlotType` but pydantic's `use_enum_values=True`
    # stores the raw string value at runtime; `PlotType(...)` normalizes either form.
    plot_type = PlotType(plot_type)
    if plot_type not in _BUILDERS:
        raise NotImplementedError(f"ECharts backend does not support {plot_type} yet")
    return _BUILDERS[plot_type]


def serialize(plot: "Plot[Any, Any]") -> Dict[str, Any]:
    """
    `{"renderer": "svg"|"canvas", "datasets": [{"layout": <option skeleton>}, ...]}`,
    parallel to `plot.model_dump()["datasets"]`. The skeleton has no `series`, no axis
    `data` arrays, and no formatter functions, so it must always be plain-JSON-safe.

    Plot types not yet ported to ECharts return `{"unsupported": <plot type>}` instead
    of raising, so a single unported plot doesn't crash the whole report render; the JS
    side (initPlot/renderPlot in echarts-plotting.js) renders a visible placeholder for
    these instead of a chart.
    """
    plot_type = PlotType(plot.plot_type)
    if plot_type not in _BUILDERS:
        return {"unsupported": plot_type.value}

    builder = _BUILDERS[plot_type]

    datasets = []
    max_marks = 0
    for dataset in plot.datasets:
        max_marks = max(max_marks, builder.mark_count(dataset))
        datasets.append({"layout": builder.layout_option(plot, dataset)})

    renderer = "canvas" if max_marks > config.echarts_canvas_threshold else "svg"
    result: Dict[str, Any] = {"renderer": renderer, "datasets": datasets}

    # Contract: the skeleton must never contain functions (those cross the Python/JS
    # bridge only via the `__FN__` sentinel in the SSR path, not here).
    json.dumps(result)

    return result


def get_option(plot: "Plot[Any, Any]", ds_idx: int, is_log: bool, is_pct: bool) -> Dict[str, Any]:
    """
    Full ECharts option for one dataset: skeleton + axis `data` + `series`, with no
    toolbox state (this is the static/SSR-export and test path, analogous to
    `Plot.get_figure`, not the interactive-render path, which JS rebuilds itself).
    """
    dataset = plot.datasets[ds_idx]
    builder = _builder(plot.plot_type)

    option = builder.layout_option(plot, dataset)
    # `axis_data` takes `pconfig` too (not just `dataset`) because the heatmap builder
    # needs `pconfig.cluster_switch_clustered_active` to pick the clustered vs.
    # unclustered category lists, and can return MULTIPLE `(axis_key, data)` entries
    # (heatmaps have both x and y category axes; bar has only y, line/scatter have none).
    axis_entries = builder.axis_data(dataset, plot.pconfig)
    if axis_entries is not None:
        for axis_key, data in axis_entries:
            option[axis_key]["data"] = data
    option["series"] = builder.series(dataset, plot.pconfig, is_pct)

    for axis_name in plot.axis_controlled_by_switches:
        echarts_axis_name = "xAxis" if axis_name == "xaxis" else "yAxis"
        axis_option = option[echarts_axis_name]
        if is_log:
            axis_option["type"] = "log"
        if is_pct:
            pct_range = dataset.pct_range.get(axis_name, {})
            axis_option["min"] = pct_range.get("min", 0)
            axis_option["max"] = pct_range.get("max", 100)

    return option
