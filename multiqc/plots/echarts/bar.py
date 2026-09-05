"""
ECharts option builder for `PlotType.BAR` (`multiqc/plots/bargraph.py`).

MultiQC bar plots are horizontal: samples on the (inverted) category axis (y),
category values on the value axis (x). See `templates/default/src/js/plots/bar.js`
for the Plotly-JS equivalent this mirrors.
"""

from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from multiqc.plots.bargraph import BarPlotConfig, Dataset
from multiqc.plots.echarts.converter import bands_and_lines, convert_layout, trailing_bands_lines_series, zoom_option

if TYPE_CHECKING:
    from multiqc.plots.plot import Plot


def _group_axis_labels(dataset: Dataset) -> List[str]:
    """
    Unique `dataset.group_labels` values, in order of first appearance. `group_labels`
    (from `pconfig.sample_groups`, see `bargraph.py::_reorder_by_groups`) is one entry
    PER ROW of `dataset.samples`/`cat.data` (a sample can repeat once per group it
    belongs to, e.g. riboWaltz's per-region bars); these unique values become the
    yAxis categories for grouped bars, one per group instead of one per row.
    """
    assert dataset.group_labels is not None
    seen: Dict[str, None] = {}
    for label in dataset.group_labels:
        seen.setdefault(label, None)
    return list(seen.keys())


def _grouped_series(dataset: Dataset, is_pct: bool) -> List[Dict[str, Any]]:
    """
    Grouped ("sample groups" / multicategory) bars. Plotly has a native multicategory
    axis and groups bars via a shared `offsetgroup` per sample (see the default
    template's `plots/bar.js`, multicategory branch); ECharts has no equivalent axis
    type, so the pragmatic port is: one ECharts `stack` id per sample (from
    `dataset.offset_groups`, falling back to the sample name), positioned at its
    `group_labels` slot on the category axis. Categories within one sample's stack sum
    together (matching how Plotly's shared-offsetgroup traces are drawn on top of each
    other); different samples (different stack ids) are placed side by side by ECharts
    automatically, which is what makes the groups visually distinct.
    """
    assert dataset.group_labels is not None
    group_labels = dataset.group_labels
    offset_groups = dataset.offset_groups or {}
    samples = dataset.samples

    group_axis = _group_axis_labels(dataset)
    group_index = {g: i for i, g in enumerate(group_axis)}

    rows_by_sample: Dict[str, List[int]] = {}
    for row_idx, sample_name in enumerate(samples):
        rows_by_sample.setdefault(sample_name, []).append(row_idx)

    result: List[Dict[str, Any]] = []
    for sample_name, row_indices in rows_by_sample.items():
        stack_id = offset_groups.get(sample_name, sample_name)
        for cat in dataset.cats:
            values = cat.data_pct if is_pct else cat.data
            data: List[Optional[float]] = [None] * len(group_axis)
            for row_idx in row_indices:
                data[group_index[group_labels[row_idx]]] = values[row_idx]
            result.append(
                {
                    "type": "bar",
                    "name": cat.name,
                    "stack": stack_id,
                    "data": data,
                    "itemStyle": {"color": cat.color},
                    "barCategoryGap": "30%",
                }
            )
    return result


def layout_option(plot: "Plot[Any, Any]", dataset: Dataset) -> Dict[str, Any]:
    """
    Full ECharts option skeleton for one bar-plot dataset, minus `series` and axis
    `data` arrays (samples are toolbox-dependent; the JS renderer fills `yAxis.data`).
    """
    layout_ir = plot.layout_ir
    option = convert_layout(layout_ir, dataset.layout)

    option["yAxis"]["type"] = "category"
    # Plotly's default (non-inverted) category axis puts the first sample at the
    # bottom and the last at the top; `dataset.samples` is already pre-reversed by
    # `bargraph.py::Dataset.create` for that convention ("Need to reverse samples as
    # the bar plot will show them reversed"). ECharts' default category axis follows
    # the same bottom-to-top convention, so leaving `inverse` unset (false) here
    # reproduces Plotly's top-to-bottom sample order exactly; `inverse: True` was
    # flipping it (verified against a live FastQC "Sequence Counts" render in both
    # engines with agent-browser).
    option["tooltip"]["trigger"] = "axis"
    # Plotly-style click+drag box-zoom (POLISH.md #17), X (value) axis only: the y-axis
    # is a sample-name category list already managed by MultiQC's own sidebar toolbox
    # (hide/highlight), see `zoom_option`.
    option.update(zoom_option(x=True, y=False))

    # JS reads this to decide whether series stack (default) or sit side by side.
    option["_mqc"] = {"barmode": layout_ir.barmode}

    # `y_bands`/`y_lines` would target the category (sample) axis here, which is
    # meaningless for a numeric threshold, so they are dropped (see
    # `converter.bands_and_lines` docstring for the empirical Plotly evidence).
    bands_lines = bands_and_lines(plot.pconfig, dataset.dconfig, include_y=False)
    if bands_lines:
        option["_mqc"]["bandsLines"] = bands_lines

    return option


def series(dataset: Dataset, pconfig: BarPlotConfig, is_pct: bool) -> List[Dict[str, Any]]:
    """
    One `{"type": "bar"}` series per category, or (grouped bars) per sample x category,
    plus (if configured) a trailing silent series carrying the static bands/lines
    markArea/markLine from `trailing_bands_lines_series`. Mirrors the stacking logic in
    `BarPlot.initialize` (`multiqc/plots/bargraph.py`).
    """
    if dataset.group_labels:
        result = _grouped_series(dataset, is_pct)
    else:
        barmode = pconfig.stacking
        if barmode is None:  # legacy: non-default None means "group"
            barmode = "group"
        if barmode == "normal":  # legacy alias
            barmode = "relative"
        is_group = barmode == "group"

        result = []
        for cat in dataset.cats:
            data = cat.data_pct if is_pct else cat.data
            s: Dict[str, Any] = {
                "type": "bar",
                "name": cat.name,
                "data": list(data),
                "itemStyle": {"color": cat.color},
                "barCategoryGap": "30%",
            }
            if not is_group:
                s["stack"] = "total"
            result.append(s)

    # See the `include_y=False` comment in `layout_option`: bar's category axis makes
    # y_bands/y_lines meaningless, so they are dropped here too.
    trailing = trailing_bands_lines_series(pconfig, dataset.dconfig, include_y=False)
    if trailing:
        result.append(trailing)

    return result


def axis_data(dataset: Dataset, pconfig: BarPlotConfig) -> List[Tuple[str, List[str]]]:
    """
    `[("yAxis", categories)]`: bar plots are always horizontal, samples on `yAxis`
    (or, for grouped bars, the unique group labels instead of raw sample rows).
    """
    if dataset.group_labels:
        return [("yAxis", _group_axis_labels(dataset))]
    return [("yAxis", list(dataset.samples))]


def mark_count(dataset: Dataset) -> int:
    return len(dataset.samples) * len(dataset.cats)
