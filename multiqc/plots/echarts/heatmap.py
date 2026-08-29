"""
ECharts option builder for `PlotType.HEATMAP` (`multiqc/plots/heatmap.py`).

Unlike bar/line/scatter, a heatmap dataset can have BOTH axes carry toolbox-dependent
sample names (`pconfig.xcats_samples`/`ycats_samples`) and a "clustered" variant of every
axis/data field (`dataset.rows_clustered`/`xcats_clustered`/`ycats_clustered`, toggled by
`pconfig.cluster_switch_clustered_active`); see
`templates/echarts/src/js/plots/heatmap.js` for the interactive (toolbox-aware) JS
counterpart of `series()`/`axis_data()` below.

Colorscale conversion (documented risk, `multiqc-echarts-exploration/BUILD_PLAN.md`
section 4 "Colorscale conversion"): Plotly's `colorscale` is a list of
`(stop_position_0_to_1, color)` pairs, interpolated at arbitrary stop positions. ECharts'
`visualMap.inRange.color` only accepts a flat color list, linearly interpolated at EVENLY
spaced stops across [min, max]. We drop the stop positions and keep only the ordered color
list. For the default 11-stop MultiQC colorscale (evenly spaced at 0.0, 0.1, ..., 1.0) this
is an exact match; for a user-supplied `pconfig.colstops` with uneven stop spacing, the
ECharts rendering will be close but not pixel-identical to Plotly. Accepted per the build
plan ("get sign-off on 'close enough'").
"""

from typing import Any, Dict, List, Tuple

from multiqc.plots.echarts.converter import convert_layout
from multiqc.plots.heatmap import Dataset, HeatmapConfig, HeatmapPlot


def _active_rows_and_cats(dataset: Dataset, cluster_active: bool) -> Tuple[List[List[Any]], List[str], List[str]]:
    """
    Resolve the clustered vs. unclustered rows/xcats/ycats for the current cluster-switch
    state, mirroring `EchartsHeatmapPlot.prepData`'s clustered-variant selection (minus the
    toolbox hide/rename filtering, which only applies in the interactive JS path).
    """
    if cluster_active and dataset.rows_clustered is not None:
        rows = dataset.rows_clustered
        xcats = list(dataset.xcats_clustered) if dataset.xcats_clustered is not None else list(dataset.xcats)
        ycats = list(dataset.ycats_clustered) if dataset.ycats_clustered is not None else list(dataset.ycats)
    else:
        rows = dataset.rows
        xcats = list(dataset.xcats)
        ycats = list(dataset.ycats)
    return rows, xcats, ycats


def _colorscale_colors(dataset: Dataset) -> List[str]:
    colorscale = dataset.trace_params.get("colorscale") or []
    colors = [str(color) for _, color in colorscale]
    if dataset.trace_params.get("reversescale"):
        colors = list(reversed(colors))
    return colors


def layout_option(plot: "HeatmapPlot", dataset: Dataset) -> Dict[str, Any]:
    """
    Full ECharts option skeleton for one heatmap dataset, minus `series` and axis `data`
    arrays (both x and y categories are potentially sample names, toolbox-dependent; the
    JS renderer fills `xAxis.data`/`yAxis.data`, see `EchartsHeatmapPlot.buildSeries`/
    `applyOptionOverrides`).
    """
    option = convert_layout(plot.layout, dataset.layout)

    option["xAxis"]["type"] = "category"
    option["yAxis"]["type"] = "category"
    # First sample/row at the top, matching Plotly's `yaxis.autorange = "reversed"`.
    option["yAxis"]["inverse"] = True
    option["tooltip"]["trigger"] = "item"

    zmin = dataset.trace_params.get("zmin")
    zmax = dataset.trace_params.get("zmax")
    vmin = zmin if zmin is not None else plot.min
    vmax = zmax if zmax is not None else plot.max
    if vmin is None or vmax is None:
        # Defensive fallback only: every non-empty heatmap has at least one numeric cell,
        # so `plot.min`/`plot.max` are set in practice (see `HeatmapPlot.create`).
        vmin, vmax = 0.0, 1.0

    option["visualMap"] = {
        "min": vmin,
        "max": vmax,
        "calculable": True,
        "orient": "horizontal",
        "left": "center",
        "bottom": 0,
        "show": bool(dataset.trace_params.get("showscale", True)),
        "inRange": {"color": _colorscale_colors(dataset)},
    }

    if plot.square:
        option["_mqc"] = {"square": True}

    return option


def series(dataset: Dataset, pconfig: HeatmapConfig, is_pct: bool) -> List[Dict[str, Any]]:
    """
    One `{"type": "heatmap"}` series with `data: [[xi, yi, val], ...]`. This is the
    SSR/get_option (non-toolbox) path; the interactive path is
    `EchartsHeatmapPlot.buildSeries()` (`templates/echarts/src/js/plots/heatmap.js`).

    `is_pct` is accepted for dispatch-signature parity with `bar.series`; heatmaps never
    enable the percentage switch, so it is unused here.
    """
    rows, xcats, _ycats = _active_rows_and_cats(dataset, pconfig.cluster_switch_clustered_active)

    # Mirrors `HeatmapPlot.create`'s own decision to enable Plotly `texttemplate` cell
    # labels: reuse that decision (rather than re-deriving the cell-count threshold here)
    # so the two backends never disagree on when labels are shown.
    show_labels = "texttemplate" in dataset.trace_params
    decimals = pconfig.tt_decimals

    data: List[Any] = []
    for yi, row in enumerate(rows):
        for xi in range(len(xcats)):
            val = row[xi] if xi < len(row) else None
            if val is None:
                continue
            if show_labels and isinstance(val, (int, float)):
                data.append(
                    {
                        "value": [xi, yi, val],
                        "label": {"show": True, "formatter": f"{val:.{decimals}f}"},
                    }
                )
            else:
                data.append([xi, yi, val])

    return [{"type": "heatmap", "data": data}]


def axis_data(dataset: Dataset, pconfig: HeatmapConfig) -> List[Tuple[str, List[str]]]:
    """`[("xAxis", xcats), ("yAxis", ycats)]` for the current cluster-switch state."""
    _, xcats, ycats = _active_rows_and_cats(dataset, pconfig.cluster_switch_clustered_active)
    return [("xAxis", xcats), ("yAxis", ycats)]


def mark_count(dataset: Dataset) -> int:
    return len(dataset.ycats) * len(dataset.xcats)
