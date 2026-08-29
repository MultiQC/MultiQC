"""
Shared Plotly-layout -> ECharts-option-skeleton conversion.

`convert_layout` reads the existing backend-agnostic Plotly `go.Layout` metadata
(titles, ranges, log/pct state) that every plot already builds, and produces the
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

import plotly.graph_objects as go  # type: ignore

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


def _axis_name(axis: Any) -> Optional[str]:
    if axis.title is not None and axis.title.text:
        return str(axis.title.text)
    return None


def _axis_type(axis: Any) -> str:
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


def _convert_axis(axis: Any, *, scale: bool = False) -> Dict[str, Any]:
    """
    `scale=True` is how a caller asks for Plotly-style autorange on a value axis: by
    default ECharts widens a value axis to always include 0, where Plotly's line/
    scatter/box traces auto-range tightly to the data extent instead (bar traces keep
    the 0 baseline, see `bar.py`, which never passes `scale=True`). Setting `scale: true`
    tells ECharts to fit the axis to the data instead of forcing in 0, matching Plotly.
    Only meaningful for a numeric axis, and only when nothing already pins the minimum
    (an explicit Plotly `range`/`autorangeoptions` min is honored as-is below).
    """
    axis_option: Dict[str, Any] = {
        "type": _axis_type(axis),
        "name": _axis_name(axis),
        "nameLocation": "middle",
        "nameGap": 30,
    }

    minval: Optional[float] = None
    maxval: Optional[float] = None
    if axis.range:
        minval, maxval = axis.range[0], axis.range[1]
    elif axis.autorangeoptions is not None:
        minval = axis.autorangeoptions.minallowed
        maxval = axis.autorangeoptions.maxallowed
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


def convert_layout(
    layout: go.Layout,
    dataset_layout: Dict[str, Any],
    *,
    scale_x: bool = False,
    scale_y: bool = False,
) -> Dict[str, Any]:
    """
    Merge `dataset_layout` (a per-dataset Plotly layout fragment, `BaseDataset.layout`)
    onto a copy of `layout` (the shared `Plot.layout`), mirroring `Plot.get_figure`'s
    `layout.update(**dataset.layout)` (`multiqc/plots/plot.py`), then convert the result
    into the shared part of an ECharts option skeleton.

    `scale_x`/`scale_y` are forwarded to `_convert_axis` as its `scale` flag, letting a
    per-type builder opt a value axis into Plotly-style data-fitted autorange instead of
    ECharts' default of always including 0 (see `_convert_axis`). Bar plots leave both
    False (bars anchor at 0); line/scatter/box pass True for their value axis/axes.
    """
    if layout is None:
        raise ValueError("converter.convert_layout: layout must not be None")

    merged = go.Layout(layout.to_plotly_json())
    merged.update(**dataset_layout)

    title_text = merged.title.text if merged.title is not None else None
    show_legend = bool(merged.showlegend)

    return {
        "animation": False,
        "title": {**_convert_title(title_text), "left": "center"},
        # ECharts' own grid defaults (left/right: "10%", top/bottom: 60) are sized for
        # a chart with no `containLabel` compensation; stacking `containLabel` on top of
        # those defaults double-counts the label space and produces much bigger margins
        # than Plotly's automargin (which sizes margins to the label content only, plus
        # a small gap). Small fixed insets here let `containLabel` do the same job:
        # grow only as far as the actual tick/axis-name text needs. The legend (to the
        # right of the grid, see `legend` below) isn't covered by `containLabel` at all,
        # so `right` needs a bit more room to keep it clear of the plot when shown, like
        # Plotly's default template (which also puts the legend on the right).
        "grid": {"left": 8, "right": 160 if show_legend else 16, "top": 64, "bottom": 8, "containLabel": True},
        "xAxis": _convert_axis(merged.xaxis, scale=scale_x),
        "yAxis": _convert_axis(merged.yaxis, scale=scale_y),
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
        "tooltip": {"confine": True},
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
