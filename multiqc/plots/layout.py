"""
Neutral, backend-agnostic layout IR.

This is the thin slice of plot-layout metadata that every rendering backend
branches on: title, per-axis type/range/log/ticksuffix/title, height/width,
legend visibility, and bar mode. It deliberately carries NOTHING Plotly- or
ECharts-specific (no colorway, hoverlabel, spikes, tickfont, margins, gridcolor,
templates): each backend owns its own rendering config and reconstructs it from
this IR plus the raw dataset fields.

Base plot modules build a `LayoutIR` (instead of a Plotly `go.Layout`); the
Plotly backend (`multiqc/plots/plotly/`) turns it into a `go.Layout`, and the
ECharts backend (`multiqc/plots/echarts/`) reads it directly. Bands/lines stay on
`pconfig` (both backends read them from there), so they are not part of the IR.
"""

from typing import Any, Dict, Literal, Optional, Tuple

from pydantic import BaseModel

AxisType = Literal["linear", "log", "category"]


class AxisIR(BaseModel):
    """One axis, in neutral terms. Empty `title`/`range` means "unset" (backend default)."""

    title: Optional[str] = None
    type: AxisType = "linear"
    # Explicit [min, max]. Either bound may be None to leave that end auto-ranged.
    range: Optional[Tuple[Optional[float], Optional[float]]] = None
    # Plotly-style soft autorange bounds (`autorangeoptions.minallowed/maxallowed`):
    # the data is fitted but never shrinks past these. Distinct from `range`, which
    # pins the view exactly.
    minallowed: Optional[float] = None
    maxallowed: Optional[float] = None
    ticksuffix: str = ""


class LayoutIR(BaseModel):
    """The neutral layout for one plot (shared across its datasets)."""

    title: Optional[str] = None  # may hold HTML like <br>/<sup>; split downstream by each backend
    height: Optional[int] = None
    width: Optional[int] = None
    showlegend: bool = False
    barmode: Optional[str] = None  # "stack" | "group" | "relative" | "overlay"
    xaxis: AxisIR = AxisIR()
    yaxis: AxisIR = AxisIR()

    def merged_with(self, override: "LayoutIR") -> "LayoutIR":
        """
        Return a copy of `self` with the explicitly-set fields of `override` applied on
        top, recursing into the axes. Mirrors the per-dataset `go.Layout.update(**fragment)`
        merge the Plotly path does today (`Plot.get_figure`): only keys the caller actually
        set on `override` win; unset fields keep `self`'s value.
        """
        merged = self.model_dump()
        _deep_update(merged, override.model_dump(exclude_unset=True))
        return LayoutIR(**merged)


def _deep_update(base: Dict[str, Any], override: Dict[str, Any]) -> None:
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _deep_update(base[key], value)
        else:
            base[key] = value
