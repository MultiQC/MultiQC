"""
ECharts SSR theme for static (flat/exported) plots.

Ports the light theme values from `multiqc-echarts-exploration/validation/echarts-themes.js`
to Python, so both the JS bundle (browser) and the SSR bootstrap in `static_export.py`
render from the same palette. Static export never switches theme (flat plots don't have a
dark-mode toggle), so only the light theme is needed here.

`FONT_FAMILY` is a SINGLE family name, not a comma-separated CSS stack: ECharts writes
`font-family` straight into the SVG `style` attribute without escaping inner quotes, so a
stack like `Inter, "Helvetica Neue", Arial` produces malformed SVG that resvg refuses to
parse (`expected a whitespace not 'H'`). See `validation/RESULTS.md` section 4.
"""

from typing import Any, Dict, List

FONT_FAMILY = "DejaVu Sans"

_PALETTE: List[str] = [
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
]

_TEXT_COLOR = "#333333"
_AXIS_COLOR = "#cccccc"
_SPLIT_COLOR = "#eeeeee"


def _axis_cfg() -> Dict[str, Any]:
    return {
        "axisLine": {"lineStyle": {"color": _AXIS_COLOR}},
        "axisTick": {"lineStyle": {"color": _AXIS_COLOR}},
        "axisLabel": {"color": _TEXT_COLOR},
        "splitLine": {"lineStyle": {"color": [_SPLIT_COLOR]}},
    }


LIGHT_THEME: Dict[str, Any] = {
    "color": _PALETTE,
    "backgroundColor": "#ffffff",
    "textStyle": {"fontFamily": FONT_FAMILY},
    "title": {"textStyle": {"color": _TEXT_COLOR}},
    "legend": {"textStyle": {"color": _TEXT_COLOR}},
    "tooltip": {"textStyle": {"fontFamily": FONT_FAMILY}},
    "categoryAxis": _axis_cfg(),
    "valueAxis": _axis_cfg(),
    "logAxis": _axis_cfg(),
    "timeAxis": _axis_cfg(),
}
