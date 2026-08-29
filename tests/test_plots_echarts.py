"""
Serialization contract tests for the ECharts backend (`multiqc/plots/echarts/`).

Uses the same bar plot as `tests/test_plots.py::test_barplot`. The autouse `reset`
fixture in `tests/conftest.py` resets `config`/`report` after every test.
"""

import math
from typing import Any, Dict

import multiqc
from multiqc import config, report
from multiqc.core.update_config import ClConfig
from multiqc.plots import bargraph, box, echarts, heatmap, linegraph, scatter, table, violin
from multiqc.plots.table_object import ColumnDict
from multiqc.types import Anchor, PlotType


def _make_bar_plot():
    return bargraph.plot(
        {
            "Sample0": {},
            "Sample1": {"Cat1": 1},
            "Sample2": {"Cat1": 1, "Cat2": 1},
            "Sample3": {"Cat1": 1, "Cat2": 1, "Cat3": 1},
        },
        ["Cat1", "Cat2"],
        bargraph.BarPlotConfig(id="bargraph", title="Test: Bar Graph"),
    )


def _make_grouped_bar_plot():
    # Mirrors the real-world usage in multiqc/modules/ribowaltz/ribowaltz.py: one
    # sample_groups group per region, each containing one row per sample; rows for the
    # same sample share an offset_group id ("sample1"/"sample2") so they align across
    # regions.
    data = {
        "sample1_5utr": {"Frame 0": 10, "Frame 1": 5, "Frame 2": 2},
        "sample2_5utr": {"Frame 0": 8, "Frame 1": 4, "Frame 2": 1},
        "sample1_cds": {"Frame 0": 20, "Frame 1": 3, "Frame 2": 1},
        "sample2_cds": {"Frame 0": 18, "Frame 1": 2, "Frame 2": 1},
    }
    sample_groups = {
        "5utr": [["sample1_5utr", "sample1"], ["sample2_5utr", "sample2"]],
        "cds": [["sample1_cds", "sample1"], ["sample2_cds", "sample2"]],
    }
    return bargraph.plot(
        data,
        ["Frame 0", "Frame 1", "Frame 2"],
        bargraph.BarPlotConfig(
            id="bargraph_grouped",
            title="Test: Grouped Bar Graph",
            sample_groups=sample_groups,
        ),
    )


def _make_line_plot():
    return linegraph.plot(
        {
            "Sample1": {0: 1.0, 1: 2.0},
            "Sample2": {0: 3.0, 1: 4.0, 2: 5.0},
        },
        linegraph.LinePlotConfig(id="linegraph", title="Test: Line Graph"),
    )


def _make_categorical_line_plot():
    return linegraph.plot(
        {
            "Sample1": {"a": 1.0, "b": 2.0},
            "Sample2": {"a": 3.0, "b": 4.0},
        },
        linegraph.LinePlotConfig(id="linegraph_categorical", title="Test: Categorical Line Graph", categories=True),
    )


def _make_line_plot_with_bands():
    return linegraph.plot(
        {"Sample1": {0: 1.0, 1: 2.0}},
        linegraph.LinePlotConfig(
            id="linegraph_bands",
            title="Test: Line Graph With Bands",
            y_bands=[{"from": 0, "to": 1, "color": "#009500", "opacity": 0.13}],
        ),
    )


def _make_scatter_plot():
    return scatter.plot(
        {
            "Sample1": [{"x": 1, "y": 2, "color": "rgb(255,0,0)", "marker_size": 12}],
            "Sample2": [{"x": 3, "y": 4}],
        },
        scatter.ScatterConfig(id="scatterplot", title="Test: Scatter"),
    )


def _make_bar_plot_with_bands():
    return bargraph.plot(
        {"Sample1": {"Cat1": 1}, "Sample2": {"Cat1": 2}},
        ["Cat1"],
        bargraph.BarPlotConfig(
            id="bargraph_bands",
            title="Test: Bar Graph With Bands",
            x_bands=[{"from": 0, "to": 1, "color": "#009500", "opacity": 0.13}],
            x_lines=[{"value": 1.5, "label": "threshold"}],
            # Meaningless for bar's category (sample) yAxis: must be dropped, not crash.
            y_bands=[{"from": 0, "to": 1, "color": "#ff0000"}],
            y_lines=[{"value": 0.5, "label": "y_ignored"}],
        ),
    )


def _make_scatter_plot_with_bands():
    return scatter.plot(
        {"Sample1": [{"x": 1, "y": 2}], "Sample2": [{"x": 3, "y": 4}]},
        scatter.ScatterConfig(
            id="scatterplot_bands",
            title="Test: Scatter With Bands",
            x_bands=[{"from": 0, "to": 1, "color": "#009500", "opacity": 0.13}],
            y_lines=[{"value": 2.5, "label": "y_threshold"}],
        ),
    )


def _make_large_scatter_plot(n_points: int):
    return scatter.plot(
        {f"Sample{i}": [{"x": i, "y": i}] for i in range(n_points)},
        scatter.ScatterConfig(id="scatterplot_large", title="Test: Large Scatter"),
    )


def _make_heatmap_plot(display_values=False):
    return heatmap.plot(
        data=[[1, 2], [3, 4]],
        xcats=["Cat1", "Cat2"],
        ycats=["Sample1", "Sample2"],
        pconfig=heatmap.HeatmapConfig(id="echarts_heatmap", title="Test: Heatmap", display_values=display_values),
    )


def _make_clustered_heatmap_plot():
    # Same fixture as test_plots.py::test_heatmap_clustering_produces_reordered_data:
    # designed so clustering actually reorders rows/cats (a 2x2 or uniform matrix may not).
    data = [
        [10, 0, 11],
        [0, 10, 0],
        [11, 0, 10],
    ]
    return heatmap.plot(
        data=data,
        xcats=["X1", "X2", "X3"],
        ycats=["A", "B", "C"],
        pconfig=heatmap.HeatmapConfig(
            id="echarts_heatmap_clustered",
            title="Test: Clustered Heatmap",
            cluster_rows=True,
            cluster_cols=True,
        ),
    )


def _make_box_plot():
    return box.plot(
        {"Sample1": [1.0, 2.0, 3.0, 4.0], "Sample2": [5.0, 6.0, 7.0, 8.0]},
        box.BoxPlotConfig(id="boxplot", title="Test: Box"),
    )


def _make_stats_box_plot():
    return box.plot(
        {
            "Sample1": {"min": 1.0, "q1": 2.0, "median": 3.0, "q3": 4.0, "max": 5.0},
            "Sample2": {"min": 2.0, "q1": 3.0, "median": 4.0, "q3": 5.0, "max": 6.0},
        },
        box.BoxPlotConfig(id="boxplot_stats", title="Test: Box Stats"),
    )


def _make_violin_plot():
    data = {f"Sample{i}": {"metric_a": float(i), "metric_b": float(i * 2), "metric_c": float(i % 5)} for i in range(10)}
    headers: Dict[str, ColumnDict] = {
        "metric_a": {"title": "Metric A"},
        "metric_b": {"title": "Metric B"},
        "metric_c": {"title": "Metric C"},
    }
    return violin.plot(data, headers, table.TableConfig(id="violinplot", title="Test: Violin"))


def test_interactive_plot_adds_echarts_key_when_engine_is_echarts():
    config.plotting_engine = "echarts"
    plot = _make_bar_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" in dumped
    assert dumped["echarts"]["renderer"] == "svg"

    skeleton = dumped["echarts"]["datasets"][0]["layout"]
    assert skeleton["animation"] is False
    assert skeleton["yAxis"]["type"] == "category"
    assert "series" not in skeleton


def test_interactive_plot_title_strips_html():
    config.plotting_engine = "echarts"
    plot = bargraph.plot(
        {
            "Sample1": {"Cat1": 1},
        },
        ["Cat1"],
        bargraph.BarPlotConfig(id="bargraph_title", title="Foo Bar<br><sup>42 things</sup>"),
    )
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    skeleton = dumped["echarts"]["datasets"][0]["layout"]
    assert skeleton["title"]["text"] == "Foo Bar"
    assert skeleton["title"]["subtext"] == "42 things"


def test_interactive_plot_unsupported_plot_type_does_not_raise(monkeypatch):
    # All PlotTypes with a figure are ported now (bar/line/scatter/heatmap/box/violin);
    # exercise the generic "unsupported" fallback path directly instead, by pretending
    # box plots aren't ported, so a single unported type still can't crash the report.
    monkeypatch.delitem(echarts._BUILDERS, PlotType.BOX)
    config.plotting_engine = "echarts"
    plot = _make_box_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert dumped["echarts"] == {"unsupported": "box plot"}


def test_interactive_plot_adds_echarts_key_for_line_plot():
    config.plotting_engine = "echarts"
    plot = _make_line_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" in dumped

    skeleton = dumped["echarts"]["datasets"][0]["layout"]
    assert skeleton["animation"] is False
    assert skeleton["xAxis"]["type"] == "value"
    # No slider (matches Plotly, which has no mini-plot strip); "inside" stays for
    # drag-pan/pinch-zoom, with mouse-wheel zoom disabled so it doesn't hijack page
    # scrolling. See `multiqc/plots/echarts/line.py::layout_option`.
    assert skeleton["dataZoom"] == [{"type": "inside", "zoomOnMouseWheel": False, "moveOnMouseWheel": False}]
    assert "series" not in skeleton


def test_interactive_plot_omits_echarts_key_by_default():
    assert config.plotting_engine == "plotly"
    plot = _make_bar_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" not in dumped


def test_get_option_series_length_matches_categories():
    plot = _make_bar_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    assert len(option["series"]) == len(dataset.cats)
    for cat, series_option in zip(dataset.cats, option["series"]):
        assert series_option["data"] == cat.data


def test_get_option_is_pct_uses_data_pct():
    plot = _make_bar_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=True)

    dataset = plot.datasets[0]
    for cat, series_option in zip(dataset.cats, option["series"]):
        # NaN != NaN, so compare positionally while tolerating NaN entries.
        for actual, expected in zip(series_option["data"], cat.data_pct):
            if math.isnan(expected):
                assert math.isnan(actual)
            else:
                assert actual == expected


def test_get_option_grouped_bar_produces_stacked_series_per_sample():
    # FIX-NEEDED #3 in multiqc-echarts-exploration/PARITY.md: grouped bars
    # (`dataset.group_labels` set, e.g. by `pconfig.sample_groups`) used to throw.
    plot = _make_grouped_bar_plot()
    dataset = plot.datasets[0]
    assert dataset.group_labels is not None

    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    unique_groups = list(dict.fromkeys(dataset.group_labels))
    assert len(unique_groups) == 2  # "5utr" and "cds"
    assert option["yAxis"]["data"] == unique_groups

    # One series per (unique sample, category): ECharts has no native multicategory
    # axis, so each sample becomes a `stack` id positioned at its group's slot instead
    # of Plotly's shared `offsetgroup`.
    unique_samples = list(dict.fromkeys(dataset.samples))
    assert len(option["series"]) == len(unique_samples) * len(dataset.cats)

    stack_ids = {s["stack"] for s in option["series"]}
    assert dataset.offset_groups is not None
    assert stack_ids == set(dataset.offset_groups.values())

    for s in option["series"]:
        assert len(s["data"]) == len(unique_groups)
        # Every row has a value in exactly one group slot for this dataset (each
        # sample appears once per region); the rest are None placeholders.
        assert sum(1 for v in s["data"] if v is not None) == 1


def test_get_option_is_log_sets_switch_controlled_axis_type():
    plot = _make_bar_plot()
    assert plot.axis_controlled_by_switches == ["xaxis"]

    option = echarts.get_option(plot, ds_idx=0, is_log=True, is_pct=False)
    assert option["xAxis"]["type"] == "log"


def test_get_option_bar_value_axis_is_not_scaled():
    # Bars must anchor at 0 (matches Plotly, which pins an explicit `minallowed=0` for
    # bar plots, see `bargraph.py`); `scale: true` is for line/scatter/box only.
    plot = _make_bar_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert option["xAxis"]["min"] == 0
    assert "scale" not in option["xAxis"]


def test_get_option_line_series_matches_lines():
    plot = _make_line_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    assert len(option["series"]) == len(dataset.lines)
    for line, series_option in zip(dataset.lines, option["series"]):
        assert series_option["type"] == "line"
        assert series_option["name"] == line.name
        assert series_option["data"] == [list(p) for p in line.pairs]
        assert series_option["smooth"] is False


def test_get_option_line_value_axis_has_no_static_data():
    plot = _make_line_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert "data" not in option["xAxis"]


def test_get_option_line_categorical_axis_uses_plain_values():
    plot = _make_categorical_line_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    assert option["xAxis"]["type"] == "category"
    assert option["xAxis"]["data"] == ["a", "b"]

    dataset = plot.datasets[0]
    for line, series_option in zip(dataset.lines, option["series"]):
        assert series_option["data"] == [y for _, y in line.pairs]


def test_get_option_line_value_axes_are_scaled():
    # Non-categorical x-axis and y-axis are both value axes for a line plot: `scale:
    # true` fits them to the data extent instead of ECharts' default of forcing 0 in,
    # matching Plotly's autorange (POLISH.md #5).
    plot = _make_line_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert option["xAxis"]["scale"] is True
    assert option["yAxis"]["scale"] is True


def test_get_option_line_categorical_axis_is_not_scaled():
    # A categorical x-axis has no meaningful "scale to data" behavior; only the
    # (still value-typed) y-axis gets it.
    plot = _make_categorical_line_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert "scale" not in option["xAxis"]
    assert option["yAxis"]["scale"] is True


def test_get_option_line_is_log_sets_yaxis():
    plot = _make_line_plot()
    assert plot.axis_controlled_by_switches == ["yaxis"]

    option = echarts.get_option(plot, ds_idx=0, is_log=True, is_pct=False)
    assert option["yAxis"]["type"] == "log"


def test_get_option_line_y_bands_produce_silent_markarea_series():
    plot = _make_line_plot_with_bands()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    marker_series = [s for s in option["series"] if s.get("silent")]
    assert len(marker_series) == 1
    assert marker_series[0]["markArea"]["data"] == [
        [{"yAxis": 0, "itemStyle": {"color": "#009500", "opacity": 0.13}}, {"yAxis": 1}]
    ]

    # The skeleton (interactive JS path) carries the same payload under `_mqc.bandsLines`.
    skeleton = echarts.serialize(plot)["datasets"][0]["layout"]
    assert skeleton["_mqc"]["bandsLines"]["markArea"] == marker_series[0]["markArea"]


def test_get_option_bar_x_bands_and_lines_produce_silent_marker_series():
    plot = _make_bar_plot_with_bands()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    marker_series = [s for s in option["series"] if s.get("silent")]
    assert len(marker_series) == 1
    # x_bands/x_lines land on the value axis (ECharts "xAxis"): bar is horizontal, so
    # xAxis is the value axis (see converter.bands_and_lines docstring).
    assert marker_series[0]["markArea"]["data"] == [
        [{"xAxis": 0, "itemStyle": {"color": "#009500", "opacity": 0.13}}, {"xAxis": 1}]
    ]
    assert marker_series[0]["markLine"]["data"][0]["xAxis"] == 1.5

    # y_bands/y_lines target bar's category (sample) axis, which is meaningless for a
    # numeric threshold: dropped cleanly, not crashed and not emitted.
    assert "yAxis" not in marker_series[0]["markArea"]["data"][0][0]
    assert all("yAxis" not in entry for entry in marker_series[0]["markLine"]["data"])

    skeleton = echarts.serialize(plot)["datasets"][0]["layout"]
    assert skeleton["_mqc"]["bandsLines"]["markArea"] == marker_series[0]["markArea"]


def test_get_option_scatter_x_bands_and_y_lines_produce_silent_marker_series():
    plot = _make_scatter_plot_with_bands()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    marker_series = [s for s in option["series"] if s.get("silent")]
    assert len(marker_series) == 1
    assert marker_series[0]["markArea"]["data"] == [
        [{"xAxis": 0, "itemStyle": {"color": "#009500", "opacity": 0.13}}, {"xAxis": 1}]
    ]
    assert marker_series[0]["markLine"]["data"] == [
        {
            "yAxis": 2.5,
            "lineStyle": {"color": None, "width": 2, "type": "solid"},
            "label": {"formatter": "y_threshold"},
        }
    ]

    skeleton = echarts.serialize(plot)["datasets"][0]["layout"]
    assert skeleton["_mqc"]["bandsLines"]["markLine"] == marker_series[0]["markLine"]


def test_interactive_plot_adds_echarts_key_for_scatter_plot():
    config.plotting_engine = "echarts"
    plot = _make_scatter_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" in dumped
    assert dumped["echarts"]["renderer"] == "svg"

    skeleton = dumped["echarts"]["datasets"][0]["layout"]
    assert skeleton["animation"] is False
    assert skeleton["tooltip"]["trigger"] == "item"
    assert "series" not in skeleton


def test_get_option_scatter_value_axes_are_scaled():
    # Scatter has no category axis: both x and y are value axes, so both get `scale:
    # true` (Plotly-style data-fitted autorange, not forced-0), matching line/box.
    plot = _make_scatter_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert option["xAxis"]["scale"] is True
    assert option["yAxis"]["scale"] is True


def test_get_option_scatter_series_has_one_series_with_value_items():
    plot = _make_scatter_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    assert len(option["series"]) == 1
    series_option = option["series"][0]
    assert series_option["type"] == "scatter"
    assert len(series_option["data"]) == len(dataset.points)

    for point, item in zip(dataset.points, series_option["data"]):
        assert item["value"] == [point["x"], point["y"]]
        assert item["name"] == point["name"]

    # Sample1's explicit color/marker_size are carried through per-item.
    sample1_item = series_option["data"][0]
    assert sample1_item["symbolSize"] == 12
    assert sample1_item["itemStyle"]["color"] == "rgb(255,0,0)"


def test_get_option_scatter_axis_has_no_static_data():
    plot = _make_scatter_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert "data" not in option["xAxis"]
    assert "data" not in option["yAxis"]


def test_scatter_mark_count_is_point_count():
    plot = _make_scatter_plot()
    assert echarts.scatter.mark_count(plot.datasets[0]) == 2


def test_serialize_scatter_uses_svg_renderer_below_threshold():
    config.plotting_engine = "echarts"
    config.echarts_canvas_threshold = 3000
    plot = _make_large_scatter_plot(10)
    result = echarts.serialize(plot)
    assert result["renderer"] == "svg"


def test_serialize_scatter_uses_canvas_renderer_above_threshold():
    config.plotting_engine = "echarts"
    config.echarts_canvas_threshold = 3000
    plot = _make_large_scatter_plot(3001)
    result = echarts.serialize(plot)
    assert result["renderer"] == "canvas"


def test_interactive_plot_adds_echarts_key_for_heatmap_plot():
    config.plotting_engine = "echarts"
    plot = _make_heatmap_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" in dumped
    assert dumped["echarts"]["renderer"] == "svg"

    skeleton = dumped["echarts"]["datasets"][0]["layout"]
    assert skeleton["animation"] is False
    assert skeleton["xAxis"]["type"] == "category"
    assert skeleton["yAxis"]["type"] == "category"
    assert "series" not in skeleton
    assert "data" not in skeleton["xAxis"]
    assert "data" not in skeleton["yAxis"]

    # visualMap: converted from the Plotly colorscale stop list (BUILD_PLAN.md "Colorscale
    # conversion" risk: stop positions are dropped, only the ordered color list survives).
    visual_map = skeleton["visualMap"]
    assert visual_map["calculable"] is True
    assert isinstance(visual_map["inRange"]["color"], list)
    assert len(visual_map["inRange"]["color"]) > 0
    assert visual_map["min"] == 1
    assert visual_map["max"] == 4


def test_get_option_heatmap_builds_xyz_cell_data():
    plot = _make_heatmap_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    assert option["xAxis"]["data"] == ["Cat1", "Cat2"]
    assert option["yAxis"]["data"] == ["Sample1", "Sample2"]

    series_option = option["series"][0]
    assert series_option["type"] == "heatmap"
    assert sorted(series_option["data"]) == [[0, 0, 1], [0, 1, 3], [1, 0, 2], [1, 1, 4]]


def test_heatmap_mark_count_is_rows_times_cols():
    plot = _make_heatmap_plot()
    assert echarts.heatmap.mark_count(plot.datasets[0]) == 4


def test_get_option_heatmap_cell_labels_when_display_values_enabled():
    plot = _make_heatmap_plot(display_values=True)
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    series_option = option["series"][0]
    for item in series_option["data"]:
        assert isinstance(item, dict)
        assert item["label"]["show"] is True
        _xi, _yi, val = item["value"]
        assert item["label"]["formatter"] == f"{val:.2f}"


def test_get_option_heatmap_clustered_switch_uses_clustered_categories():
    plot = _make_clustered_heatmap_plot()
    dataset = plot.datasets[0]
    assert dataset.rows_clustered is not None  # sanity: clustering actually happened

    plot.pconfig.cluster_switch_clustered_active = True
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    assert option["xAxis"]["data"] == list(dataset.xcats_clustered)
    assert option["yAxis"]["data"] == list(dataset.ycats_clustered)

    plot.pconfig.cluster_switch_clustered_active = False
    option_unclustered = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert option_unclustered["xAxis"]["data"] == list(dataset.xcats)
    assert option_unclustered["yAxis"]["data"] == list(dataset.ycats)


# GOLDEN quartile test: this fixed input + expected five-number/outlier values is the
# cross-language contract asserted here AND mirrored in a comment block at the top of
# `multiqc/templates/echarts/src/js/plots/box.js`. The JS `fiveNumberSummary()`/
# `outliers()` port must produce the same output for the same input (checked manually:
# there is no JS unit-test runner in this repo).
_GOLDEN_BOX_VALUES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100]
_GOLDEN_FIVE_NUMBER_SUMMARY = [1.0, 3.5, 6.0, 8.5, 10.0]
_GOLDEN_OUTLIERS = [100.0]


def test_box_five_number_summary_golden_values():
    assert echarts.box.five_number_summary(_GOLDEN_BOX_VALUES) == _GOLDEN_FIVE_NUMBER_SUMMARY


def test_box_outliers_golden_values():
    assert echarts.box.outliers(_GOLDEN_BOX_VALUES) == _GOLDEN_OUTLIERS


def test_interactive_plot_adds_echarts_key_for_box_plot():
    config.plotting_engine = "echarts"
    plot = _make_box_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" in dumped
    assert dumped["echarts"]["renderer"] == "svg"

    skeleton = dumped["echarts"]["datasets"][0]["layout"]
    assert skeleton["animation"] is False
    assert skeleton["yAxis"]["type"] == "category"
    assert skeleton["yAxis"]["inverse"] is True
    assert skeleton["tooltip"]["trigger"] == "item"
    assert "series" not in skeleton


def test_get_option_box_series_has_five_number_data():
    plot = _make_box_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    boxplot_series = next(s for s in option["series"] if s["type"] == "boxplot")
    assert len(boxplot_series["data"]) == len(dataset.samples)
    for values, entry in zip(dataset.data, boxplot_series["data"]):
        assert entry["value"] == echarts.box.five_number_summary(values)
        # Semi-transparent fill + solid border, so the median (drawn via borderColor)
        # is visible against the box body: POLISH.md #6.
        assert entry["itemStyle"]["color"].startswith("rgba(")
        assert entry["itemStyle"]["color"].endswith(",0.5)")
        assert entry["itemStyle"]["borderColor"].startswith("rgb(")


def test_get_option_box_outliers_become_scatter_series():
    plot = _make_box_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    scatter_series = [s for s in option["series"] if s["type"] == "scatter"]
    # Neither sample's 4 evenly-spaced values produce a Tukey outlier.
    total_outliers = sum(len(echarts.box.outliers(values)) for values in dataset.data)
    assert total_outliers == 0
    assert scatter_series == []


def test_get_option_box_axis_has_no_static_sample_data():
    plot = _make_box_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert option["yAxis"]["data"] == list(plot.datasets[0].samples)
    assert "data" not in option["xAxis"]


def test_get_option_box_value_axis_is_scaled():
    # Box plots are horizontal: xAxis is the value axis and gets `scale: true`
    # (Plotly-style data-fitted autorange, POLISH.md #7); yAxis is the (inverted)
    # category axis and is untouched by scale.
    plot = _make_box_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert option["xAxis"]["scale"] is True
    assert option["yAxis"]["type"] == "category"
    assert "scale" not in option["yAxis"]


def test_get_option_box_stats_data_uses_precomputed_values_directly():
    plot = _make_stats_box_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    assert dataset.is_stats_data is True
    boxplot_series = next(s for s in option["series"] if s["type"] == "boxplot")
    for stats, entry in zip(dataset.data, boxplot_series["data"]):
        assert entry["value"] == [stats["min"], stats["q1"], stats["median"], stats["q3"], stats["max"]]
    # No raw values to derive outliers from.
    assert not any(s["type"] == "scatter" for s in option["series"])


def test_box_mark_count_is_samples_plus_outliers():
    plot = _make_box_plot()
    dataset = plot.datasets[0]
    expected = len(dataset.samples) + sum(len(echarts.box.outliers(values)) for values in dataset.data)
    assert echarts.box.mark_count(dataset) == expected

    outlier_dataset = box.Dataset(
        **{**dataset.__dict__, "samples": ["S1"], "data": [_GOLDEN_BOX_VALUES], "is_stats_data": False}
    )
    assert echarts.box.mark_count(outlier_dataset) == 1 + len(_GOLDEN_OUTLIERS)


# GOLDEN kde() test: this fixed input + expected output is the cross-language contract
# asserted here AND mirrored in a comment block at the top of `multiqc/plots/echarts/violin.py`
# (and, in Task 2.2, at the top of `templates/echarts/src/js/plots/violin.js`). The JS
# `kde()` port must reproduce these same values for the same input.
_GOLDEN_KDE_VALUES = [1.0, 2.0, 3.0, 4.0, 5.0]
_GOLDEN_KDE_XS = [1.0, 3.0, 5.0]
_GOLDEN_KDE_DENSITIES = [0.15916497933387785, 0.1802710624663249, 0.15916497933387785]


def test_kde_golden_values():
    assert echarts.violin.kde(_GOLDEN_KDE_VALUES, _GOLDEN_KDE_XS) == _GOLDEN_KDE_DENSITIES


def _renderitem_body(series_option):
    return series_option["renderItem"]["body"]


def test_interactive_plot_adds_echarts_key_for_violin_plot():
    config.plotting_engine = "echarts"
    plot = _make_violin_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" in dumped
    assert dumped["echarts"]["renderer"] == "svg"

    dataset_entry = dumped["echarts"]["datasets"][0]
    skeleton = dataset_entry["layout"]
    assert skeleton["animation"] is False
    assert skeleton["xAxis"]["type"] == "value"
    assert skeleton["xAxis"]["show"] is False
    assert skeleton["yAxis"]["type"] == "value"
    assert skeleton["yAxis"]["inverse"] is True
    assert skeleton["yAxis"]["min"] == -0.5
    assert skeleton["yAxis"]["max"] == 2.5  # 3 metrics -> rows 0, 1, 2
    assert skeleton["yAxis"]["axisLabel"]["formatter"]["__FN__"] is True
    assert "series" not in skeleton

    # The skeleton itself must stay plain-JSON-safe (the `__FN__` sentinel is a dict, not
    # an actual function): `serialize()` already asserts this globally via `json.dumps`.

    # VIOLIN EXCEPTION: precomputed KDE polygons ride along per-dataset, keyed by metric.
    violins = dataset_entry["violins"]
    assert set(violins.keys()) == {"metric_a", "metric_b", "metric_c"}
    for payload in violins.values():
        assert len(payload["poly"]) == 2 * echarts.violin.N_BINS
        assert len(payload["range"]) == 2
        assert payload["range"][0] < payload["range"][1]


def test_get_option_violin_series_counts():
    plot = _make_violin_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    n_metrics = len(dataset.metrics)
    assert n_metrics == 3

    kde_series = [s for s in option["series"] if s["type"] == "custom" and "polygon" in _renderitem_body(s)]
    annotation_series = [s for s in option["series"] if s["type"] == "custom" and "polygon" not in _renderitem_body(s)]
    scatter_series = [s for s in option["series"] if s["type"] == "scatter"]

    assert len(kde_series) == n_metrics
    assert len(annotation_series) == n_metrics
    assert len(scatter_series) == n_metrics
    assert len(option["series"]) == 3 * n_metrics


def test_get_option_violin_scatter_items_shape():
    plot = _make_violin_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    scatter_by_name = {s["name"]: s for s in option["series"] if s["type"] == "scatter"}
    assert set(scatter_by_name) == set(dataset.metrics)

    for metric, series_option in scatter_by_name.items():
        n_values = len(dataset.violin_value_by_sample_by_metric[metric])
        assert len(series_option["data"]) == n_values
        for item in series_option["data"]:
            assert list(item.keys()) == ["value", "name"]
            normx, row_idx = item["value"]
            assert 0.0 <= normx <= 1.0
            assert isinstance(row_idx, int)
            assert item["name"] in dataset.all_samples


def test_get_option_violin_axis_has_no_static_data():
    plot = _make_violin_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)
    assert "data" not in option["xAxis"]
    assert "data" not in option["yAxis"]


def test_violin_mark_count_is_total_scatter_points():
    plot = _make_violin_plot()
    dataset = plot.datasets[0]
    expected = sum(len(dataset.violin_value_by_sample_by_metric[m]) for m in dataset.metrics)
    assert echarts.violin.mark_count(dataset) == expected == 30  # 3 metrics * 10 samples


def test_serialize_violin_includes_violins_payload():
    plot = _make_violin_plot()
    result = echarts.serialize(plot)

    violins = result["datasets"][0]["violins"]
    dataset = plot.datasets[0]
    assert set(violins.keys()) == set(dataset.metrics)
    for metric in dataset.metrics:
        payload = violins[metric]
        assert len(payload["poly"]) == 2 * echarts.violin.N_BINS
        for x, y in payload["poly"]:
            assert isinstance(x, float)
            assert isinstance(y, float)
        assert payload["range"] == list(echarts.violin._metric_range(echarts.violin._numeric_values(dataset, metric)))


############################################
# Task 3.2 (multiqc-echarts-exploration/BUILD_PLAN.md "test suite integration"):
# engine-parametrized tests, a box sort-by-median unit test, and end-to-end pipeline
# tests.


def _verify_rendered_for_engine(plot, plotting_engine: str) -> Dict[str, Any]:
    """
    Cloned subset of `tests/test_plots.py::_verify_rendered`, parametrized over both
    backends via the `plotting_engine` fixture (`tests/conftest.py`): render the plot and
    check the payload contract for whichever engine is active.
    """
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))
    dumped = report.plot_data[plot.anchor]
    if plotting_engine == "echarts":
        assert "echarts" in dumped
        payload = dumped["echarts"]
        assert "unsupported" not in payload
        assert payload["renderer"] in ("svg", "canvas")
        assert len(payload["datasets"]) == len(plot.datasets)
    else:
        assert "echarts" not in dumped
    return dumped


def test_engine_parametrized_bar(plotting_engine):
    _verify_rendered_for_engine(_make_bar_plot(), plotting_engine)


def test_engine_parametrized_line(plotting_engine):
    _verify_rendered_for_engine(_make_line_plot(), plotting_engine)


def test_engine_parametrized_scatter(plotting_engine):
    _verify_rendered_for_engine(_make_scatter_plot(), plotting_engine)


def test_engine_parametrized_heatmap(plotting_engine):
    _verify_rendered_for_engine(_make_heatmap_plot(), plotting_engine)


def test_engine_parametrized_box(plotting_engine):
    _verify_rendered_for_engine(_make_box_plot(), plotting_engine)


def test_engine_parametrized_violin(plotting_engine):
    _verify_rendered_for_engine(_make_violin_plot(), plotting_engine)


def test_box_sort_by_median_produces_median_sorted_series():
    """
    D3 (`multiqc-echarts-exploration/PARITY.md`): the sort-by-median toggle is wired up
    only on the JS side (the inherited, engine-neutral `prepData()` swaps in
    `dataset.data_sorted`/`samples_sorted` when the sort switch is active,
    `templates/echarts/src/js/plots/box.js`); no shipped module sets `sort_by_median=True`
    and nothing exercises that sorted path. `Dataset.create` (`multiqc/plots/box.py`)
    always keeps `dataset.data`/`samples` as the alphabetical order and stashes the
    ascending-by-median order separately in `data_sorted`/`samples_sorted`. Simulate
    exactly what the JS toggle does (swap those in as the active dataset fields) and
    drive the result through the real `echarts.get_option()` used by the SSR/flat export
    path, proving the ECharts box builder (which only ever reads `dataset.data`/
    `samples`) is order-agnostic and renders a genuinely median-sorted series.
    """
    data = {
        "SampleB": [1.0, 2.0, 3.0],  # median 2
        "SampleA": [10.0, 20.0, 30.0],  # median 20
        "SampleC": [100.0, 200.0, 300.0],  # median 200
    }
    plot = box.plot(
        data,
        box.BoxPlotConfig(id="boxplot_sorted", title="Test: Box Sorted", sort_by_median=True),
    )
    assert isinstance(plot, box.BoxPlot)
    dataset = plot.datasets[0]
    assert dataset.samples_sorted == ["SampleB", "SampleA", "SampleC"]
    assert dataset.data_sorted is not None

    sorted_dataset = box.Dataset(**{**dataset.__dict__, "samples": dataset.samples_sorted, "data": dataset.data_sorted})
    plot.datasets = [sorted_dataset]

    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    assert option["yAxis"]["data"] == ["SampleB", "SampleA", "SampleC"]
    boxplot_series = next(s for s in option["series"] if s["type"] == "boxplot")
    medians = [entry["value"][2] for entry in boxplot_series["data"]]
    assert medians == [2.0, 20.0, 200.0]
    assert medians == sorted(medians)


############################################
# End-to-end: the full `multiqc.run()` pipeline against a single, self-contained module
# data file (fastp's SAMPLE.json; no directory scan needed) so this stays fast. fastp
# emits bar, line, and violin plots, so this also touches most of the ported plot types.


def test_echarts_template_end_to_end_interactive(data_dir, tmp_path):
    result = multiqc.run(
        data_dir / "modules" / "fastp" / "SAMPLE.json",
        cfg=ClConfig(run_modules=["fastp"], template="echarts", output_dir=tmp_path),
        return_html=True,
    )

    assert result.sys_exit_code == 0
    assert result.html_content is not None
    # The echarts template's footer credit (multiqc/templates/echarts/footer.html):
    # proves the echarts bundle/template was actually used, not just the default one.
    assert "Apache ECharts" in result.html_content

    assert len(report.plot_data) > 0
    for dumped in report.plot_data.values():
        assert "echarts" in dumped
        assert "unsupported" not in dumped["echarts"]


def test_echarts_template_end_to_end_flat(data_dir, tmp_path):
    result = multiqc.run(
        data_dir / "modules" / "fastp" / "SAMPLE.json",
        cfg=ClConfig(run_modules=["fastp"], template="echarts", plots_force_flat=True, output_dir=tmp_path),
        return_html=True,
    )

    assert result.sys_exit_code == 0
    assert result.html_content is not None
    # Static export path (SSR + resvg PNG rasterization): flat plots are embedded as
    # base64 PNGs, mirroring the Plotly flat-plot contract.
    assert "data:image/png;base64," in result.html_content
