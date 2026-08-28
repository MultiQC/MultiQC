"""
Shared Plotly-layout -> ECharts-option-skeleton conversion.

`convert_layout` reads the existing backend-agnostic Plotly `go.Layout` metadata
(titles, ranges, log/pct state) that every plot already builds, and produces the
part of an ECharts `option` dict that is common to every plot type: `animation`,
`title`, `grid`, `xAxis`/`yAxis`, `legend`, `tooltip`. Per-type builders (e.g.
`multiqc/plots/echarts/bar.py`) call this first, then add type-specific fields.

This module never adds `series`, axis `data` arrays, or formatter functions: see
the "ECharts model->JSON contract" in `multiqc-echarts-exploration/BUILD_PLAN.md`.
"""

from typing import Any, Dict, Optional

import plotly.graph_objects as go  # type: ignore


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


def _convert_axis(axis: Any) -> Dict[str, Any]:
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
    if maxval is not None:
        axis_option["max"] = maxval

    if axis.ticksuffix:
        # String template, not a function: functions only cross the Python/JS bridge
        # in the SSR path, via the `__FN__` sentinel (see static_export.py, later task).
        axis_option["axisLabel"] = {"formatter": "{value}" + str(axis.ticksuffix)}

    return axis_option


def convert_layout(layout: go.Layout, dataset_layout: Dict[str, Any]) -> Dict[str, Any]:
    """
    Merge `dataset_layout` (a per-dataset Plotly layout fragment, `BaseDataset.layout`)
    onto a copy of `layout` (the shared `Plot.layout`), mirroring `Plot.get_figure`'s
    `layout.update(**dataset.layout)` (`multiqc/plots/plot.py`), then convert the result
    into the shared part of an ECharts option skeleton.
    """
    if layout is None:
        raise ValueError("converter.convert_layout: layout must not be None")

    merged = go.Layout(layout.to_plotly_json())
    merged.update(**dataset_layout)

    title_text = merged.title.text if merged.title is not None else None

    return {
        "animation": False,
        "title": {"text": title_text, "left": "center"},
        "grid": {"containLabel": True},
        "xAxis": _convert_axis(merged.xaxis),
        "yAxis": _convert_axis(merged.yaxis),
        "legend": {
            "show": bool(merged.showlegend),
            "orient": "horizontal",
            "bottom": 0,
            "left": "center",
        },
        "tooltip": {"confine": True},
    }
