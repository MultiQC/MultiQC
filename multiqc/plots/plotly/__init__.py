"""
Plotly rendering backend.

Turns the neutral `LayoutIR` (`multiqc/plots/layout.py`) into Plotly `go.Layout`
objects, and back. This package OWNS everything Plotly-specific (the MultiQC
template, margins, hover styling, spikes, tickfonts, the legend placement): base
plot modules no longer build `go.Layout` directly.

`ir_to_layout` / `layout_to_ir` are the two halves of the bridge. `layout_to_ir`
extracts the neutral slice from an existing `go.Layout` (used while migrating base
modules onto the IR, and to prove `ir_to_layout` reconstructs faithfully); once a
module builds the IR natively, only `ir_to_layout` runs.
"""

from multiqc.plots.plotly.layout_builder import ir_to_layout, layout_to_ir

__all__ = ["ir_to_layout", "layout_to_ir"]
