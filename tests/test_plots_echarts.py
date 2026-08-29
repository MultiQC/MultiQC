"""
Serialization contract tests for the ECharts backend (`multiqc/plots/echarts/`).

Uses the same bar plot as `tests/test_plots.py::test_barplot`. The autouse `reset`
fixture in `tests/conftest.py` resets `config`/`report` after every test.
"""

import math

from multiqc import config, report
from multiqc.plots import bargraph, box, echarts, heatmap, linegraph, scatter, table, violin
from multiqc.types import Anchor


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
    # Violin is not ported to ECharts yet (Phase 2): used to verify the "unsupported
    # plot type" fallback still works now that bar/line/scatter/heatmap/box are supported.
    return violin.plot(
        data={"Sample1": {"x": 1, "y": 2}, "Sample2": {"x": 3, "y": 4}},
        headers={"x": {"title": "Metric X"}},
        pconfig=table.TableConfig(id="violinplot", title="Test: Violin"),
    )


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


def test_interactive_plot_unsupported_plot_type_does_not_raise():
    config.plotting_engine = "echarts"
    plot = _make_violin_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert dumped["echarts"] == {"unsupported": "violin plot"}


def test_interactive_plot_adds_echarts_key_for_line_plot():
    config.plotting_engine = "echarts"
    plot = _make_line_plot()
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))

    dumped = report.plot_data[plot.anchor]
    assert "echarts" in dumped

    skeleton = dumped["echarts"]["datasets"][0]["layout"]
    assert skeleton["animation"] is False
    assert skeleton["xAxis"]["type"] == "value"
    assert skeleton["dataZoom"] == [{"type": "inside"}, {"type": "slider"}]
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


def test_get_option_is_log_sets_switch_controlled_axis_type():
    plot = _make_bar_plot()
    assert plot.axis_controlled_by_switches == ["xaxis"]

    option = echarts.get_option(plot, ds_idx=0, is_log=True, is_pct=False)
    assert option["xAxis"]["type"] == "log"


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
    for values, five_number in zip(dataset.data, boxplot_series["data"]):
        assert five_number == echarts.box.five_number_summary(values)


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


def test_get_option_box_stats_data_uses_precomputed_values_directly():
    plot = _make_stats_box_plot()
    option = echarts.get_option(plot, ds_idx=0, is_log=False, is_pct=False)

    dataset = plot.datasets[0]
    assert dataset.is_stats_data is True
    boxplot_series = next(s for s in option["series"] if s["type"] == "boxplot")
    for stats, five_number in zip(dataset.data, boxplot_series["data"]):
        assert five_number == [stats["min"], stats["q1"], stats["median"], stats["q3"], stats["max"]]
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
