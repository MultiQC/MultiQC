"""
`LayoutIR` <-> Plotly `go.Layout`.

`ir_to_layout` reconstructs the base MultiQC Plotly layout (the fixed scaffolding
that `Plot.initialize` builds today: template, title anchor, automargined axes,
margins, hover label, legend-when-flat) from the neutral IR plus constants.
`layout_to_ir` does the inverse extraction of the neutral slice. Keeping the fixed
Plotly styling here, driven off the small IR, is what lets base modules stop
importing Plotly.
"""

from typing import Any, Dict, Optional

import plotly.graph_objects as go  # type: ignore

from multiqc.plots.layout import AxisIR, AxisType, LayoutIR

# Fixed base-layout constants, lifted verbatim from `Plot.initialize` so that
# `ir_to_layout` reproduces its output byte-for-byte (guarded by tests/test_layout_ir.py).
_MARGIN = dict(pad=5, t=50, r=15, b=65, l=60)
_HOVERLABEL_NAMELENGTH = -1
_FLAT_LEGEND = dict(orientation="h", yanchor="top", y=-0.15, xanchor="center", x=0.5)


def _axis_ir_to_dict(axis: AxisIR) -> Dict[str, Any]:
    """Neutral axis -> Plotly axis dict. `type` is emitted only when it is not the Plotly
    auto default ("linear" is left unset so numeric axes keep Plotly's auto-detect, matching
    `Plot.initialize`, which only ever sets `type` to "log")."""
    d: Dict[str, Any] = {"automargin": True}
    if axis.type != "linear":
        d["type"] = axis.type
    if axis.title is not None:
        d["title"] = {"text": axis.title}
    if axis.range is not None:
        d["range"] = list(axis.range)
    if axis.minallowed is not None or axis.maxallowed is not None:
        d["autorangeoptions"] = {"minallowed": axis.minallowed, "maxallowed": axis.maxallowed}
    if axis.ticksuffix:
        d["ticksuffix"] = axis.ticksuffix
    return d


def ir_to_layout(ir: LayoutIR, *, flat: bool = False, extra: Optional[Dict[str, Any]] = None) -> go.Layout:
    """Build the base MultiQC `go.Layout` from the neutral IR. `flat` adds the horizontal
    bottom legend the static/flat path uses (interactive reports get no `legend` block).
    `extra` is a per-plot plotly-json fragment (`Plot.plotly_layout_extra`) replayed on top
    with `go.Layout.update`, reproducing the layout modules used to mutate directly."""
    # Lazy import: the MultiQC Plotly template lives in plot.py, which imports this package.
    from multiqc.plots.plot import get_multiqc_plotly_template

    kwargs: Dict[str, Any] = dict(
        template=go.layout.Template(get_multiqc_plotly_template()),
        title=go.layout.Title(text=ir.title, xanchor="center", x=0.5),
        xaxis=go.layout.XAxis(**_axis_ir_to_dict(ir.xaxis)),
        yaxis=go.layout.YAxis(**_axis_ir_to_dict(ir.yaxis)),
        height=ir.height,
        width=ir.width,
        autosize=True,
        margin=go.layout.Margin(**_MARGIN),
        hoverlabel=go.layout.Hoverlabel(namelength=_HOVERLABEL_NAMELENGTH),
        showlegend=ir.showlegend,
        legend=go.layout.Legend(**_FLAT_LEGEND) if flat else None,
    )
    if ir.barmode is not None:
        kwargs["barmode"] = ir.barmode
    layout = go.Layout(**kwargs)
    if extra:
        layout.update(**extra)
    return layout


def _axis_dict_to_ir(axis: Any) -> AxisIR:
    axis_type: AxisType = axis.type if axis.type in ("log", "category") else "linear"
    title = axis.title.text if axis.title is not None and axis.title.text else None
    range_ = tuple(axis.range) if axis.range else None  # type: ignore[arg-type]
    auto = axis.autorangeoptions
    minallowed = auto.minallowed if auto is not None else None
    maxallowed = auto.maxallowed if auto is not None else None
    return AxisIR(
        title=title,
        type=axis_type,
        range=range_,
        minallowed=minallowed,
        maxallowed=maxallowed,
        ticksuffix=str(axis.ticksuffix) if axis.ticksuffix else "",
    )


def layout_to_ir(layout: go.Layout) -> LayoutIR:
    """Extract the neutral slice from an existing Plotly `go.Layout`."""
    if layout is None:
        raise ValueError("layout_to_ir: layout must not be None")
    title: Optional[str] = layout.title.text if layout.title is not None else None
    return LayoutIR(
        title=title,
        height=layout.height,
        width=layout.width,
        showlegend=bool(layout.showlegend),
        barmode=layout.barmode,
        xaxis=_axis_dict_to_ir(layout.xaxis),
        yaxis=_axis_dict_to_ir(layout.yaxis),
    )
