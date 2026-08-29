"""
ECharts option builder for `PlotType.BAR` (`multiqc/plots/bargraph.py`).

MultiQC bar plots are horizontal: samples on the (inverted) category axis (y),
category values on the value axis (x). See `templates/default/src/js/plots/bar.js`
for the Plotly-JS equivalent this mirrors.
"""

from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from multiqc.plots.bargraph import BarPlotConfig, Dataset
from multiqc.plots.echarts.converter import convert_layout

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
    option = convert_layout(plot.layout, dataset.layout)

    option["yAxis"]["type"] = "category"
    option["yAxis"]["inverse"] = True
    option["tooltip"]["trigger"] = "axis"

    # JS reads this to decide whether series stack (default) or sit side by side.
    option["_mqc"] = {"barmode": plot.layout.barmode}

    return option


def series(dataset: Dataset, pconfig: BarPlotConfig, is_pct: bool) -> List[Dict[str, Any]]:
    """
    One `{"type": "bar"}` series per category, or (grouped bars) per sample x category.
    Mirrors the stacking logic in `BarPlot.initialize` (`multiqc/plots/bargraph.py`).
    """
    if dataset.group_labels:
        return _grouped_series(dataset, is_pct)

    barmode = pconfig.stacking
    if barmode is None:  # legacy: non-default None means "group"
        barmode = "group"
    if barmode == "normal":  # legacy alias
        barmode = "relative"
    is_group = barmode == "group"

    result: List[Dict[str, Any]] = []
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
