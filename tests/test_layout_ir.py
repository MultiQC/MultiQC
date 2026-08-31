"""
Parity + unit tests for the neutral layout IR and the Plotly bridge.

The load-bearing check is `test_base_layout_roundtrips`: it proves the Plotly
backend's `ir_to_layout` reconstructs `Plot.initialize`'s base `go.Layout`
byte-for-byte from the small IR, so base modules can be migrated onto the IR
without changing what the Plotly renderer receives.
"""

import pytest

from multiqc.plots.bargraph import BarPlotConfig
from multiqc.plots.layout import AxisIR, LayoutIR
from multiqc.plots.plot import Plot
from multiqc.plots.plotly import ir_to_layout, layout_to_ir
from multiqc.types import Anchor, PlotType


def _base_plot(**pconfig_kwargs) -> Plot:
    pconfig = BarPlotConfig(id="ir_test", title="IR Test", **pconfig_kwargs)
    return Plot.initialize(
        plot_type=PlotType.BAR,
        pconfig=pconfig,
        anchor=Anchor("ir_test_anchor"),
        n_series_per_dataset=[1],
    )


@pytest.mark.parametrize(
    "pconfig_kwargs",
    [
        {},  # plain base
        {"height": 800, "width": 400},  # explicit size
        {"ylog": True},  # log axis -> type="log" must survive the round-trip
        {"logswitch": True, "logswitch_active": True},  # switch-driven log on the y axis
    ],
)
def test_base_layout_roundtrips(pconfig_kwargs):
    base = _base_plot(**pconfig_kwargs).layout
    rebuilt = ir_to_layout(layout_to_ir(base), flat=False)
    assert rebuilt.to_plotly_json() == base.to_plotly_json()


def test_ir_captures_the_semantic_slice():
    ir = layout_to_ir(_base_plot(height=650, ylog=True).layout)
    assert ir.title == "IR Test"
    assert ir.height == 650
    assert ir.yaxis.type == "log"
    assert ir.xaxis.type == "linear"


def test_ir_to_layout_flat_adds_bottom_legend():
    interactive = ir_to_layout(LayoutIR(title="t"), flat=False)
    flat = ir_to_layout(LayoutIR(title="t"), flat=True)
    assert interactive.legend.orientation is None
    assert flat.legend.orientation == "h"


def test_ir_to_layout_carries_axis_details():
    ir = LayoutIR(
        title="t",
        barmode="stack",
        xaxis=AxisIR(title="Reads", ticksuffix="%", minallowed=0, maxallowed=100),
        yaxis=AxisIR(type="category"),
    )
    layout = ir_to_layout(ir)
    assert layout.barmode == "stack"
    assert layout.xaxis.title.text == "Reads"
    assert layout.xaxis.ticksuffix == "%"
    assert layout.xaxis.autorangeoptions.minallowed == 0
    assert layout.xaxis.autorangeoptions.maxallowed == 100
    assert layout.yaxis.type == "category"
    # A plain "linear" axis stays unset so Plotly keeps its numeric auto-detect.
    assert layout.yaxis.type != "linear"  # category was set
    assert ir_to_layout(LayoutIR()).xaxis.type is None


def test_merged_with_overrides_only_set_fields():
    base = LayoutIR(title="base", height=500, showlegend=True, xaxis=AxisIR(title="x", ticksuffix="%"))
    override = LayoutIR(height=800, xaxis=AxisIR(title="x2"))
    merged = base.merged_with(override)
    assert merged.height == 800  # overridden
    assert merged.title == "base"  # untouched
    assert merged.showlegend is True  # untouched
    assert merged.xaxis.title == "x2"  # overridden
    assert merged.xaxis.ticksuffix == "%"  # untouched within the axis
