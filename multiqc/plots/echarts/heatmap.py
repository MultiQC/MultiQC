"""
ECharts option builder for `PlotType.HEATMAP` (`multiqc/plots/heatmap.py`).

Unlike bar/line/scatter, a heatmap dataset can have BOTH axes carry toolbox-dependent
sample names (`pconfig.xcats_samples`/`ycats_samples`) and a "clustered" variant of every
axis/data field (`dataset.rows_clustered`/`xcats_clustered`/`ycats_clustered`, toggled by
`pconfig.cluster_switch_clustered_active`); see
`templates/default/src/js/plots-echarts/heatmap.js` for the interactive (toolbox-aware) JS
counterpart of `series()`/`axis_data()` below.

Colorscale conversion: Plotly's `colorscale` is a list of `(stop_position_0_to_1, color)`
pairs, interpolated at arbitrary stop positions. ECharts' `visualMap.inRange.color` only
accepts a flat color list, linearly interpolated at EVENLY spaced stops across [min, max].
To reproduce Plotly's mapping exactly (including uneven stop spacing, e.g. FastQC's
pass/warn/fail status heatmap, and hard colour steps from duplicate/adjacent stop
positions), `_colorscale_colors` below RESAMPLES the colorscale into a fine, evenly spaced
list of colours (`_RESAMPLE_STEPS`), piecewise-linearly interpolated in RGBA at each of
those even steps. Feeding that fine, even list to `visualMap.inRange.color` makes ECharts'
own even-spacing interpolation reproduce Plotly's uneven one to within the resample
resolution.
"""

import re
from typing import Any, Dict, List, Tuple

from multiqc.plots.echarts.converter import convert_layout
from multiqc.plots.heatmap import Dataset, HeatmapConfig, HeatmapPlot

# Number of evenly spaced colours ECharts' visualMap is fed: 256 resample intervals (257
# points), so binary-fraction stop positions (0.25, 0.5, 0.125, ...), the common case for
# hand-picked colorscales like FastQC's pass/warn/fail bands, land exactly on a resample
# point instead of being split across two, and any other spacing is still visually
# indistinguishable from Plotly's own continuous interpolation.
_RESAMPLE_STEPS = 257

_HEX_RE = re.compile(r"^#([0-9a-fA-F]{3,8})$")
_RGB_RE = re.compile(r"^rgba?\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*(?:,\s*([\d.]+)\s*)?\)$", re.IGNORECASE)
# ponytail: only the CSS keywords MultiQC colorscales are known to use in practice
# (`multiqc/modules/*/`.py`s `colstops` are otherwise all `#hex`/`rgba()`); extend this
# table if a colorscale ever needs an uncommon named colour.
_NAMED_COLORS: Dict[str, str] = {
    "white": "#ffffff",
    "black": "#000000",
    "red": "#ff0000",
    "green": "#008000",
    "blue": "#0000ff",
    "yellow": "#ffff00",
    "cyan": "#00ffff",
    "aqua": "#00ffff",
    "magenta": "#ff00ff",
    "fuchsia": "#ff00ff",
    "gray": "#808080",
    "grey": "#808080",
    "orange": "#ffa500",
    "purple": "#800080",
    "navy": "#000080",
    "lime": "#00ff00",
    "maroon": "#800000",
    "olive": "#808000",
    "teal": "#008080",
    "silver": "#c0c0c0",
    "transparent": "rgba(0, 0, 0, 0)",
}

RGBA = Tuple[float, float, float, float]


def _parse_color(color: str) -> RGBA:
    """Parse a `#hex` (3/4/6/8 digit), `rgb()`/`rgba()`, or CSS-named colour string into
    `(r, g, b, a)` with r/g/b in 0..255 and a in 0..1."""
    color = _NAMED_COLORS.get(color.strip().lower(), color.strip())

    m = _RGB_RE.match(color)
    if m:
        r, g, b = (float(m.group(i)) for i in (1, 2, 3))
        a = float(m.group(4)) if m.group(4) is not None else 1.0
        return r, g, b, a

    m = _HEX_RE.match(color)
    if m:
        digits = m.group(1)
        if len(digits) in (3, 4):
            digits = "".join(c * 2 for c in digits)
        if len(digits) == 6:
            digits += "ff"
        if len(digits) == 8:
            r, g, b, a = (int(digits[i : i + 2], 16) for i in (0, 2, 4, 6))
            return float(r), float(g), float(b), a / 255.0

    raise ValueError(f"Unrecognised heatmap colorscale colour: {color!r}")


def _format_rgba(rgba: RGBA) -> str:
    r, g, b, a = rgba
    return f"rgba({round(r)}, {round(g)}, {round(b)}, {round(a, 3)})"


def _interp(c0: RGBA, c1: RGBA, frac: float) -> RGBA:
    return (
        c0[0] + (c1[0] - c0[0]) * frac,
        c0[1] + (c1[1] - c0[1]) * frac,
        c0[2] + (c1[2] - c0[2]) * frac,
        c0[3] + (c1[3] - c0[3]) * frac,
    )


def _color_at(stops: List[Tuple[float, RGBA]], t: float) -> RGBA:
    """Piecewise-linear interpolation at `t`, matching Plotly's own colorscale lookup."""
    if t <= stops[0][0]:
        return stops[0][1]
    if t >= stops[-1][0]:
        return stops[-1][1]
    for (p0, c0), (p1, c1) in zip(stops, stops[1:]):
        if t > p1:
            continue
        if p1 == p0:
            # Duplicate/adjacent stop positions: a hard colour step in the source
            # colorscale. Jump straight to the colour after the step instead of
            # blending (a blend would require dividing by a zero-width span).
            return c1
        return _interp(c0, c1, (t - p0) / (p1 - p0))
    return stops[-1][1]  # pragma: no cover (unreachable: covered by the bounds checks above)


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
    """
    Resample `dataset.trace_params["colorscale"]` (Plotly's `[(stop_position, color), ...]`,
    positions in 0..1) into `_RESAMPLE_STEPS` evenly spaced colours, so ECharts'
    even-spacing-only `visualMap.inRange.color` reproduces Plotly's piecewise-linear,
    possibly-uneven-stop mapping. See the module docstring.
    """
    colorscale = dataset.trace_params.get("colorscale") or []
    if not colorscale:
        return []

    stops = sorted(((float(pos), _parse_color(str(color))) for pos, color in colorscale), key=lambda s: s[0])

    if dataset.trace_params.get("reversescale"):
        # Plotly's own `reversescale` mirrors stop positions (`1 - pos`), not just the
        # colour order (verified against `plotly`'s rendered output for an uneven,
        # non-symmetric colorscale): reversing only the colour order at the original
        # positions is exact solely when the stops happen to be symmetric (e.g. the
        # default evenly spaced scale), which is why that simpler approach used to work.
        stops = sorted(((1.0 - pos, color) for pos, color in stops), key=lambda s: s[0])

    lo, hi = stops[0][0], stops[-1][0]
    span = hi - lo or 1.0
    return [_format_rgba(_color_at(stops, lo + span * i / (_RESAMPLE_STEPS - 1))) for i in range(_RESAMPLE_STEPS)]


def layout_option(plot: "HeatmapPlot", dataset: Dataset) -> Dict[str, Any]:
    """
    Full ECharts option skeleton for one heatmap dataset, minus `series` and axis `data`
    arrays (both x and y categories are potentially sample names, toolbox-dependent; the
    JS renderer fills `xAxis.data`/`yAxis.data`, see `EchartsHeatmapPlot.buildSeries`/
    `applyOptionOverrides`).
    """
    option = convert_layout(plot.layout_ir, dataset.layout)

    option["xAxis"]["type"] = "category"
    option["yAxis"]["type"] = "category"
    # Angle the x-axis sample labels like Plotly's own `xaxis.tickangle = 45`
    # (`HeatmapPlot.create`, `multiqc/plots/heatmap.py`), so many more (often all) of
    # them fit before ECharts' auto label-skipping kicks in, instead of most being
    # dropped as horizontal text. `grid.containLabel` (set by `convert_layout` above)
    # already grows the plot area to fit the taller rotated bounding box, so no extra
    # bottom margin is needed here.
    option["xAxis"]["axisLabel"] = {"rotate": 45}
    # First sample/row at the top, matching Plotly's `yaxis.autorange = "reversed"`.
    option["yAxis"]["inverse"] = True
    option["tooltip"]["trigger"] = "item"
    # No click+drag box-zoom (POLISH.md #17): unlike bar/box/line/scatter, ECharts'
    # toolbox dataZoomSelect brush never reliably activates for a `heatmap` series in
    # this build (confirmed empirically with agent-browser across canvas/SVG renderer,
    # square/non-square layout, and multiple activation strategies: the toolbox model
    # itself reports the cursor mode as active, but the brush's drag-to-select never
    # engages, so no range ever changes). Rather than ship a control that silently does
    # nothing, heatmap keeps its pre-existing no-zoom state; see README.md's "Known
    # gaps" section. `zoom_option` is intentionally not called here.

    zmin = dataset.trace_params.get("zmin")
    zmax = dataset.trace_params.get("zmax")
    vmin = zmin if zmin is not None else plot.min
    vmax = zmax if zmax is not None else plot.max
    if vmin is None or vmax is None:
        # Defensive fallback only: every non-empty heatmap has at least one numeric cell,
        # so `plot.min`/`plot.max` are set in practice (see `HeatmapPlot.create`).
        vmin, vmax = 0.0, 1.0

    show_scale = bool(dataset.trace_params.get("showscale", True))
    # Vertical colorbar on the right, matching Plotly's default `colorbar` placement.
    # `itemHeight` is a raw pixel count (ECharts' visualMap ignores percent strings here,
    # see `VisualMapModel.resetExtent`), so it's derived from `plot.layout.height` (the
    # same total figure height Plotly uses for its own colorbar length) minus the grid's
    # top margin and room below the bar for its own min-value label plus the x-axis tick
    # labels beneath the grid. Clamped so tiny plots don't get a negative/zero bar.
    plot_height = plot.layout_ir.height or 400
    item_height = max(100, plot_height - 64 - 90)
    option["visualMap"] = {
        "min": vmin,
        "max": vmax,
        "calculable": True,
        "orient": "vertical",
        "right": 8,
        "top": 64,
        "itemHeight": item_height,
        "show": show_scale,
        "inRange": {"color": _colorscale_colors(dataset)},
    }
    if show_scale:
        # The colorbar sits to the right of the grid at a fixed `right: 8`; it isn't
        # accounted for by `containLabel`, so widen the grid's right margin here too, or
        # the colorbar overlaps the last data column. 90 clears the colorbar's own
        # rendered width (itemWidth plus its min/max value labels, ~70px measured
        # empirically) plus a small gap.
        option["grid"]["right"] = 90

    if plot.square:
        option["_mqc"] = {"square": True}

    return option


def series(dataset: Dataset, pconfig: HeatmapConfig, is_pct: bool) -> List[Dict[str, Any]]:
    """
    One `{"type": "heatmap"}` series with `data: [[xi, yi, val], ...]`. This is the
    SSR/get_option (non-toolbox) path; the interactive path is
    `EchartsHeatmapPlot.buildSeries()` (`templates/default/src/js/plots-echarts/heatmap.js`).

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
