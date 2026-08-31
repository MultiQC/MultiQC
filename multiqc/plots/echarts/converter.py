"""
Shared neutral-layout -> ECharts-option-skeleton conversion.

`convert_layout` reads the backend-agnostic `LayoutIR` (`multiqc/plots/layout.py`)
metadata (titles, ranges, log/pct state) that every plot builds, and produces the
part of an ECharts `option` dict that is common to every plot type: `animation`,
`title`, `grid`, `xAxis`/`yAxis`, `legend`, `tooltip`. Per-type builders (e.g.
`multiqc/plots/echarts/bar.py`) call this first, then add type-specific fields.

This module never adds `series`, axis `data` arrays, or live formatter functions: see
the "ECharts model->JSON contract" in `multiqc-echarts-exploration/BUILD_PLAN.md`. It
does emit `__FN__` sentinels (plain, JSON-safe dicts) for value/log axis tick labels,
same pattern violin.py uses for its yAxis; see `_si_axis_formatter_body`.

`bands_and_lines`/`echarts_dash` also live here (not in `line.py`) since bar.py and
scatter.py need the same `x_bands`/`y_bands`/`x_lines`/`y_lines` -> markArea/markLine
conversion; line.py was the original (and still the reference) caller.
"""

import html
import json
import re
from typing import TYPE_CHECKING, Any, Dict, List, Optional

from multiqc.plots.layout import AxisIR, LayoutIR

if TYPE_CHECKING:
    from multiqc.plots.plot import PConfig

_BR_RE = re.compile(r"<br\s*/?>", re.IGNORECASE)
_TAG_RE = re.compile(r"<[^>]+>")
_WHITESPACE_RE = re.compile(r"\s+")


def _clean_title_segment(segment: str) -> str:
    unescaped = html.unescape(segment)
    stripped = _TAG_RE.sub("", unescaped)
    return _WHITESPACE_RE.sub(" ", stripped).strip()


def _convert_title(title_text: Optional[str]) -> Dict[str, Optional[str]]:
    """
    Convert a Plotly title string (which may contain HTML like `<br>` and `<sup>`)
    into ECharts `title.text`/`title.subtext` fields, which are plain text.
    """
    if not title_text:
        return {"text": None, "subtext": None}

    segments = _BR_RE.split(title_text)
    text = _clean_title_segment(segments[0])
    subtext = _clean_title_segment(" ".join(segments[1:])) if len(segments) > 1 else None
    return {"text": text, "subtext": subtext or None}


def _axis_name(axis: AxisIR) -> Optional[str]:
    return axis.title


def _axis_type(axis: AxisIR) -> str:
    if axis.type == "log":
        return "log"
    if axis.type == "category":
        return "category"
    return "value"


def _si_axis_formatter_body(suffix: str) -> str:
    """
    `__FN__` sentinel body (see module docstring / `static_export.py`) for a value/log
    axis tick label: SI-abbreviates large numbers to ~3 significant figures, trimming
    trailing zeros (450000000 -> "450M", 1500000 -> "1.5M", 12000 -> "12k"), matching
    Plotly's default axis tick style (POLISH.md #12), then appends `suffix` verbatim
    (e.g. a coverage "x" unit).

    GOLDEN CROSS-LANGUAGE CONTRACT: this algorithm must stay identical to
    `formatAxisNumber()` in `multiqc/templates/echarts/src/js/echarts-plotting.js` (same
    duplication pattern as `violin.py`'s `kde()`/JS `kde()` pair) since only the SSR path
    (`static_export.py`) ever executes this body; the interactive JS renderer always
    overwrites this sentinel with its own live `formatAxisNumber` call before
    `setOption()`, so a sentinel object never reaches ECharts in the browser.
    """
    return (
        "if (typeof v !== 'number' || !isFinite(v)) return String(v);"
        "var sign = v < 0 ? '-' : '';"
        "var abs = Math.abs(v);"
        "var units = [[1e12, 'T'], [1e9, 'G'], [1e6, 'M'], [1e3, 'k']];"
        "for (var i = 0; i < units.length; i++) {"
        "if (abs >= units[i][0]) {"
        "return sign + Number((abs / units[i][0]).toPrecision(3)) + units[i][1] + " + json.dumps(suffix) + ";"
        "}"
        "}"
        "return sign + Number(abs.toPrecision(3)) + " + json.dumps(suffix) + ";"
    )


# Gap (px) ECharts leaves between an axis's tick labels and its `name` (title), when
# the axis has one. Reused below (`_grid_inset`) so the grid's left/bottom inset can
# reserve exactly this much extra room, rather than the two numbers drifting apart.
_AXIS_NAME_GAP = 30


def _convert_axis(axis: AxisIR, *, scale: bool = False) -> Dict[str, Any]:
    """
    `scale=True` is how a caller asks for Plotly-style autorange on a value axis: by
    default ECharts widens a value axis to always include 0, where Plotly's line/
    scatter/box traces auto-range tightly to the data extent instead (bar traces keep
    the 0 baseline, see `bar.py`, which never passes `scale=True`). Setting `scale: true`
    tells ECharts to fit the axis to the data instead of forcing in 0, matching Plotly.
    Only meaningful for a numeric axis, and only when nothing already pins the minimum
    (an explicit `range`/`minallowed` on the IR is honored as-is below).
    """
    axis_option: Dict[str, Any] = {
        "type": _axis_type(axis),
        "name": _axis_name(axis),
        "nameLocation": "middle",
        "nameGap": _AXIS_NAME_GAP,
    }

    minval: Optional[float] = None
    maxval: Optional[float] = None
    if axis.range is not None:
        minval, maxval = axis.range[0], axis.range[1]
    else:
        minval = axis.minallowed
        maxval = axis.maxallowed
    if minval is not None:
        axis_option["min"] = minval
    elif scale and axis_option["type"] != "category":
        axis_option["scale"] = True
    if maxval is not None:
        axis_option["max"] = maxval

    if axis_option["type"] in ("value", "log"):
        # SI-abbreviate large tick labels (POLISH.md #12), matching Plotly's default axis
        # style. A category axis (sample names) is never numeric, so it's left untouched.
        suffix = str(axis.ticksuffix) if axis.ticksuffix else ""
        axis_option["axisLabel"] = {"formatter": {"__FN__": True, "body": _si_axis_formatter_body(suffix)}}

    return axis_option


def _grid_inset(axis_option: Dict[str, Any]) -> int:
    """
    Px inset `convert_layout` reserves on an axis's outer grid edge (`left` for yAxis,
    `bottom` for xAxis).

    `containLabel` (set on `grid` below) already grows the plot area to fit tick labels,
    but NOT the axis `name`/title (POLISH.md #14, a known ECharts limitation: `name` is
    drawn `nameGap` further out than the tick labels, and `containLabel`'s reserved space
    stops at the labels). Left at the same tiny inset used when there's no name, a
    present name gets pushed past the container edge and disappears entirely (confirmed
    on RiboTish's "Read Length Distribution" and somalier's "Relatedness" scatter, both
    losing their axis titles outright, not just clipping them short). So: reserve
    `_AXIS_NAME_GAP` (the same gap `nameGap` uses) plus one line of name text when a name
    is set, and fall back to the pre-existing tiny inset otherwise, keeping bar/box/violin/
    heatmap axes (none of which set a `name`) exactly as tight as before.
    """
    if axis_option.get("name"):
        return _AXIS_NAME_GAP + 14
    return 8


def convert_layout(
    layout: LayoutIR,
    dataset_layout: Dict[str, Any],
    *,
    scale_x: bool = False,
    scale_y: bool = False,
) -> Dict[str, Any]:
    """
    Merge `dataset_layout` (a per-dataset layout fragment, `BaseDataset.layout`, a plain
    dict) onto `layout` (the shared `Plot.layout_ir`, a neutral `LayoutIR`), mirroring
    `Plot.get_figure`'s `layout.update(**dataset.layout)` (`multiqc/plots/plot.py`), then
    convert the result into the shared part of an ECharts option skeleton. No Plotly here:
    the IR carries exactly the metadata (titles, ranges, log/pct state) ECharts needs.

    `scale_x`/`scale_y` are forwarded to `_convert_axis` as its `scale` flag, letting a
    per-type builder opt a value axis into Plotly-style data-fitted autorange instead of
    ECharts' default of always including 0 (see `_convert_axis`). Bar plots leave both
    False (bars anchor at 0); line/scatter/box pass True for their value axis/axes.
    """
    if layout is None:
        raise ValueError("converter.convert_layout: layout must not be None")

    merged = layout.merged_with(LayoutIR.from_dataset_layout(dataset_layout))

    title_text = merged.title
    show_legend = merged.showlegend

    x_axis_option = _convert_axis(merged.xaxis, scale=scale_x)
    y_axis_option = _convert_axis(merged.yaxis, scale=scale_y)

    return {
        "animation": False,
        "title": {**_convert_title(title_text), "left": "center"},
        # ECharts' own grid defaults (left/right: "10%", top/bottom: 60) are sized for
        # a chart with no `containLabel` compensation; stacking `containLabel` on top of
        # those defaults double-counts the label space and produces much bigger margins
        # than Plotly's automargin (which sizes margins to the label content only, plus
        # a small gap). Small fixed insets here let `containLabel` do the same job:
        # grow only as far as the actual tick/axis-name text needs, except `left`/`bottom`
        # widen further when the respective axis has a `name` set, so that name is never
        # cropped (see `_grid_inset`, POLISH.md #14). The legend (to the right of the
        # grid, see `legend` below) isn't covered by `containLabel` at all, so `right`
        # needs a bit more room to keep it clear of the plot when shown, like Plotly's
        # default template (which also puts the legend on the right).
        "grid": {
            "left": _grid_inset(y_axis_option),
            "right": 160 if show_legend else 16,
            "top": 64,
            "bottom": _grid_inset(x_axis_option),
            "containLabel": True,
        },
        "xAxis": x_axis_option,
        "yAxis": y_axis_option,
        "legend": {
            "show": bool(merged.showlegend),
            "type": "scroll",
            "orient": "vertical",
            "right": 8,
            "top": 56,
            "bottom": 8,
            # `right`'s 160px budget covers the icon plus a normal-length category name
            # (Plotly instead grows the margin to fit, which we can't do without a text
            # measurer on the Python side); truncate the rare much-longer label instead
            # of letting it overflow into the plot area.
            "textStyle": {"width": 110, "overflow": "truncate"},
        },
        "tooltip": {"confine": True, "axisPointer": {"type": "cross"}},
        # Plotly-style crosshair guide lines (POLISH.md #13): dashed lines on both axes
        # that track the cursor, with small axis-value label boxes. `tooltip.axisPointer`
        # above is enough for an axis-trigger tooltip (bar/line, `option["tooltip"]["trigger"]
        # = "axis"`, set per-type); an item-trigger tooltip (scatter/box/heatmap/violin) never
        # shows an axis pointer on its own, so this top-level `axisPointer` component is what
        # actually draws the cross for those types. Harmless to set for every type: it's
        # ignored by the SSR/flat PNG path (no hover there) and, for axis-trigger charts,
        # simply provides the same cross style tooltip.axisPointer already asks for.
        # Color/label theming is layered on afterwards in JS (buildCurrentOption's theme
        # step), same pattern as tooltip/axis/legend colors below it; violin.py turns this
        # off (see its layout_option) since its axes are normalized/hidden and a value
        # label here would be actively misleading.
        "axisPointer": {"type": "cross", "show": True, "triggerOn": "mousemove"},
    }


# `Series.dash` is already normalized to Plotly dash names by `convert_dash_style`
# (multiqc/plots/plot.py); ECharts `lineStyle.type` only knows "solid"/"dashed"/"dotted".
_DASH_MAP = {
    "solid": "solid",
    "dash": "dashed",
    "longdash": "dashed",
    "dot": "dotted",
    "dashdot": "dashed",
    "longdashdot": "dashed",
}


def echarts_dash(dash: Optional[str]) -> str:
    if dash is None:
        return "solid"
    return _DASH_MAP.get(dash, "solid")


def bands_and_lines(
    pconfig: "PConfig",
    *,
    include_x: bool = True,
    include_y: bool = True,
) -> Dict[str, Any]:
    """
    Static markArea/markLine payload for `pconfig.x_bands`/`y_bands`/`x_lines`/`y_lines`
    (e.g. the fastqc per-base-quality green/yellow/red zones), mirroring the Plotly
    `shapes` built by `Plot._set_x_bands_and_range`/`_set_y_bands_and_range`
    (`multiqc/plots/plot.py`). These are static (independent of toolbox state), so this
    is computed once and reused by both `layout_option` (stashed in the skeleton under
    `_mqc.bandsLines` for the interactive JS path) and `series` (the SSR/get_option path).

    `include_x`/`include_y` let a caller drop a dimension that has no meaningful target
    axis. Bar plots are horizontal (samples on a category `yAxis`, values on `xAxis`,
    see `multiqc/plots/echarts/bar.py`), so `y_bands`/`y_lines` would place a numeric
    band/line against a category axis (meaningless, would misplace or crash); bar.py
    calls this with `include_y=False`. Verified empirically against the Plotly
    reference (`BarPlot.get_figure`): a bar plot's `x_bands` land on the final Plotly
    x-axis (the value axis after bar's own xaxis/yaxis title swap), matching
    `xAxis` here; `y_bands` land on the category axis and render uselessly there.
    Line and scatter plots use both axes as meaningful value axes, so both default to
    `True`.

    Per-dataset `data_labels` overrides of bands/lines are not supported yet (parity gap,
    same class of gap as bar's sample-groups gap); only the top-level `pconfig` fields
    are read.
    """
    mark_area_data: List[List[Dict[str, Any]]] = []
    if include_y:
        for band in pconfig.y_bands or []:
            mark_area_data.append(
                [
                    {"yAxis": band.from_, "itemStyle": {"color": band.color, "opacity": band.opacity}},
                    {"yAxis": band.to},
                ]
            )
    if include_x:
        for band in pconfig.x_bands or []:
            mark_area_data.append(
                [
                    {"xAxis": band.from_, "itemStyle": {"color": band.color, "opacity": band.opacity}},
                    {"xAxis": band.to},
                ]
            )

    mark_line_data: List[Dict[str, Any]] = []
    if include_y:
        for line in pconfig.y_lines or []:
            mark_line_data.append(
                {
                    "yAxis": line.value,
                    "lineStyle": {"color": line.color, "width": line.width, "type": echarts_dash(line.dash)},
                    "label": {"formatter": line.label or ""},
                }
            )
    if include_x:
        for line in pconfig.x_lines or []:
            mark_line_data.append(
                {
                    "xAxis": line.value,
                    "lineStyle": {"color": line.color, "width": line.width, "type": echarts_dash(line.dash)},
                    "label": {"formatter": line.label or ""},
                }
            )

    result: Dict[str, Any] = {}
    if mark_area_data:
        result["markArea"] = {"silent": True, "data": mark_area_data}
    if mark_line_data:
        result["markLine"] = {"silent": True, "symbol": "none", "data": mark_line_data}
    return result


def zoom_option(*, x: bool = True, y: bool = True) -> Dict[str, Any]:
    """
    Plotly-style click+drag box-zoom (POLISH.md #17). ECharts' rubber-band box-zoom lives
    behind the toolbox `dataZoom` feature; the interactive JS renderer
    (`templates/echarts/src/js/echarts-plotting.js`) activates its cursor mode
    automatically after every render via `takeGlobalCursor`, so the user never has to
    click a toolbar button first, matching Plotly's default click-drag-to-zoom. The
    toolbox icon row itself is pushed off-canvas (`top: "150%"`): MultiQC already has its
    own sidebar toolbox, so a second visible one would just be clutter, but the feature
    still has to be declared here for `takeGlobalCursor` to find it. An `inside` dataZoom
    on the same axes holds the current zoom range (so pinch-zoom and the box-zoom above
    stay in sync); wheel zoom/pan stay off so a long report's mouse wheel keeps scrolling
    the page instead of zooming a chart under the cursor.

    `x`/`y` restrict which axis is zoomable, matching which axis is semantically
    continuous for a given plot type (see each `layout_option`'s call site: `bar.py`/
    `box.py` pass `y=False`, since their category y-axis is a sample list already managed
    by MultiQC's own sidebar toolbox (hide/highlight), not something worth a second,
    duplicate "zoom" interaction; `heatmap.py` passes both, since it has no numeric axis
    at all, both sample-category axes ARE the zoomable space; `line.py`/`scatter.py` pass
    both, since both axes carry real numeric/positional data either way).
    `violin.py` never calls this: its axes are a normalized 0..1-per-row x-trick and a
    row-index y-trick (see its module docstring), not real data coordinates, so zoom is
    disabled outright, the same reasoning POLISH.md #13 already used to disable its
    axisPointer.
    """
    x_idx = [0] if x else []
    y_idx = [0] if y else []
    # One `inside` dataZoom PER zoomable axis, not one combined entry declaring both
    # xAxisIndex and yAxisIndex together: empirically (agent-browser drag test on a
    # scatter plot), a single dataZoom spanning both axes never actually updates its
    # start/end when the toolbox's box-select zoom fires, even though the identical
    # setup with one axis per component (bar/box's x-only case) works correctly. Two
    # independent components is also what ECharts' own box-zoom examples use.
    inside: List[Dict[str, Any]] = []
    if x:
        inside.append({"type": "inside", "xAxisIndex": [0], "zoomOnMouseWheel": False, "moveOnMouseWheel": False})
    if y:
        inside.append({"type": "inside", "yAxisIndex": [0], "zoomOnMouseWheel": False, "moveOnMouseWheel": False})
    return {
        "toolbox": {
            "show": True,
            "top": "150%",
            "feature": {"dataZoom": {"show": True, "xAxisIndex": x_idx, "yAxisIndex": y_idx}},
        },
        "dataZoom": inside,
    }


def trailing_bands_lines_series(
    pconfig: "PConfig",
    *,
    include_x: bool = True,
    include_y: bool = True,
) -> Optional[Dict[str, Any]]:
    """
    `bands_and_lines` wrapped as a trailing, invisible `{"type": "line"}` series (or
    `None` if `pconfig` has no bands/lines): appended to `series()` in line.py/bar.py/
    scatter.py so it shows up in the SSR/get_option (non-toolbox) path. The interactive
    JS path carries the same `bands_and_lines` payload separately, via `layout_option`
    stashing it under `option["_mqc"]["bandsLines"]` for `Plot.bandsLinesSeries()`
    (`templates/echarts/src/js/echarts-plotting.js`) to read.
    """
    bands_lines = bands_and_lines(pconfig, include_x=include_x, include_y=include_y)
    if not bands_lines:
        return None
    return {
        "type": "line",
        "name": "",
        "data": [],
        "silent": True,
        "showSymbol": False,
        "tooltip": {"show": False},
        **bands_lines,
    }
