import sys
import tempfile
from typing import Dict, List, Union
from unittest.mock import patch

import pytest

from multiqc import config, report
from multiqc.core.exceptions import RunError
from multiqc.plots import bargraph, box, heatmap, linegraph, scatter, table, violin
from multiqc.plots.linegraph import LinePlotConfig, Series
from multiqc.plots.plot import Plot, process_batch_exports
from multiqc.plots.table_object import ColumnDict
from multiqc.types import Anchor
from multiqc.validation import ModuleConfigValidationError


@pytest.fixture(autouse=True)
def reset_config():
    """Reset config state during tests that modify global config."""
    original_boxplot_boxpoints = config.boxplot_boxpoints
    original_box_min_threshold_no_points = config.box_min_threshold_no_points
    original_box_min_threshold_outliers = config.box_min_threshold_outliers
    original_development = config.development
    original_export_plots = config.export_plots
    original_export_plot_formats = getattr(config, "export_plot_formats", None)
    original_strict = config.strict
    yield
    config.boxplot_boxpoints = original_boxplot_boxpoints
    config.box_min_threshold_no_points = original_box_min_threshold_no_points
    config.box_min_threshold_outliers = original_box_min_threshold_outliers
    config.development = original_development
    config.export_plots = original_export_plots
    if original_export_plot_formats is not None:
        config.export_plot_formats = original_export_plot_formats
    elif hasattr(config, "export_plot_formats"):
        delattr(config, "export_plot_formats")
    config.strict = original_strict


def _verify_rendered(plot) -> Plot:
    assert isinstance(plot, Plot)
    plot.add_to_report(module_anchor=Anchor("test"), section_anchor=Anchor("test"))
    assert len(report.plot_data) == 1
    assert plot.id in report.plot_data
    return plot


############################################
# Normal plot tests first, no special cases.


def test_barplot():
    _verify_rendered(
        bargraph.plot(
            {
                "Sample0": {},
                "Sample1": {"Cat1": 1},
                "Sample2": {"Cat1": 1, "Cat2": 1},
                "Sample3": {"Cat1": 1, "Cat2": 1, "Cat3": 1},
            },
            ["Cat1", "Cat2"],
            bargraph.BarPlotConfig(id="bargraph", title="Test: Bar Graph"),
        )
    )


def test_linegraph():
    dataset = {
        "Sample0": {},
        "Sample1": {0: 1, 1: 1},
        "Sample2": {0: 1, 1: 1, 2: 1},
    }
    plot = _verify_rendered(
        linegraph.plot(
            dataset,
            linegraph.LinePlotConfig(id="linegraph", title="Test: Line Graph"),
        )
    )

    assert len(report.plot_data[plot.anchor]["datasets"][0]["lines"][0]["pairs"]) == 2
    assert len(report.plot_data[plot.anchor]["datasets"][0]["lines"][1]["pairs"]) == 3


def test_table():
    _verify_rendered(
        table.plot(
            data={
                "sample1": {"x": 1, "y": 2},
                "sample2": {"x": 3, "y": 4},
            },
            headers={"x": {"title": "Metric X"}},
            pconfig=table.TableConfig(id="table", title="Table"),
        )
    )


def test_violin():
    _verify_rendered(
        violin.plot(
            data={"sample1": {"x": 1, "y": 2}, "sample2": {"x": 3, "y": 4}},
            headers={"x": {"title": "Metric X"}},
            pconfig=table.TableConfig(id="violin", title="Violin"),
        )
    )


def test_heatmap():
    _verify_rendered(
        heatmap.plot(
            data=[[1, 2], [3, 4]],
            xcats=["Cat1", "Cat2"],
            ycats=["Sample1", "Sample2"],
            pconfig=heatmap.HeatmapConfig(id="heatmap", title="Heatmap"),
        )
    )


def test_scatter():
    _verify_rendered(
        scatter.plot(
            {"Sample1": [{"x": 1, "y": 2}]},
            scatter.ScatterConfig(id="scatter", title="Scatter", xlab="X", ylab="Y"),
        )
    )


def test_boxplot():
    _verify_rendered(
        box.plot(
            {"Sample1": [1, 2, 3], "Sample2": [4, 5, 6]},
            box.BoxPlotConfig(id="box", title="Box"),
        )
    )


@pytest.mark.parametrize(
    "boxpoints_value",
    ["outliers", "all", "suspectedoutliers", False],
)
def test_boxplot_custom_boxpoints(boxpoints_value):
    """
    Test box plot with custom boxpoints configuration using global config
    """
    config.boxplot_boxpoints = boxpoints_value

    data = {
        "Sample": [
            # Tight distribution with few outliers
            30,
            31,
            32,
            33,
            34,
            35,
            36,
            37,
            38,
            39,
            40,
            # Outliers
            20,
            50,
        ],
    }
    # Test with "all" boxpoints (show all data points)
    plot_all = _verify_rendered(
        box.plot(
            data,  # type: ignore
            box.BoxPlotConfig(id=f"box_{boxpoints_value}", title=f"Box Plot - {boxpoints_value}"),
        )
    )

    plot_data_all = report.plot_data[plot_all.anchor]
    trace_params_all = plot_data_all["datasets"][0]["trace_params"]
    assert trace_params_all["boxpoints"] == boxpoints_value


def test_boxplot_dynamic_boxpoints():
    """
    Test box plot dynamic boxpoints behavior based on sample count
    """
    # Reset config to test dynamic behavior
    config.boxplot_boxpoints = None

    # Test with few samples (should show all points)
    config.box_min_threshold_no_points = 10
    config.box_min_threshold_outliers = 5

    data_few: Dict[str, List[Union[int, float]]] = {
        "Sample1": [1.0, 2.0, 3.0, 4.0, 5.0],
        "Sample2": [2.0, 3.0, 4.0, 5.0, 6.0],
    }

    plot_few = _verify_rendered(
        box.plot(
            data_few,
            box.BoxPlotConfig(id="box_few_samples", title="Box Plot - Few Samples"),
        )
    )

    plot_data_few = report.plot_data[plot_few.anchor]
    trace_params_few = plot_data_few["datasets"][0]["trace_params"]
    assert trace_params_few["boxpoints"] == "all"

    report.reset()

    # Test with many samples (should show only outliers)
    data_many: Dict[str, List[Union[int, float]]] = {f"Sample{i}": [1.0, 2.0, 3.0, 4.0, 5.0] for i in range(10)}

    plot_many = _verify_rendered(
        box.plot(
            data_many,
            box.BoxPlotConfig(id="box_many_samples", title="Box Plot - Many Samples"),
        )
    )

    plot_data_many = report.plot_data[plot_many.anchor]
    trace_params_many = plot_data_many["datasets"][0]["trace_params"]
    assert trace_params_many["boxpoints"] == "outliers"

    report.reset()

    # Test with very many samples (should show no points)
    data_very_many: Dict[str, List[Union[int, float]]] = {f"Sample{i}": [1.0, 2.0, 3.0, 4.0, 5.0] for i in range(15)}

    plot_very_many = _verify_rendered(
        box.plot(
            data_very_many,
            box.BoxPlotConfig(id="box_very_many_samples", title="Box Plot - Very Many Samples"),
        )
    )

    plot_data_very_many = report.plot_data[plot_very_many.anchor]
    trace_params_very_many = plot_data_very_many["datasets"][0]["trace_params"]
    assert trace_params_very_many["boxpoints"] is False


############################################
# Plot special cases.


def test_bar_plot_no_matching_cats():
    """
    None of the cats are matching thd data, so shouldn't produce a plot
    """
    plot_id = "test_bar_plot_no_matching_cats"

    plot = bargraph.plot(
        {"Sample1": {"Cat0": 1, "Cat1": 1}},
        ["Cat2", "Cat3"],
        {"id": plot_id, "title": "Test: Bar Graph"},
    )
    # Will return a warning message html instead of a plot:
    assert plot is None


def test_bar_plot_cats_dicts():
    """
    Advanced cats spec - dict with cat properties instead of a simple list
    """
    plot = _verify_rendered(
        bargraph.plot(
            {"Sample1": {"Cat1": 1}},
            {"Cat1": {"name": "My category"}},
            {"id": "test_bar_plot_cats_dicts", "title": "Test: Bar Graph"},
        )
    )
    assert report.plot_data[plot.anchor]["datasets"][0]["cats"][0]["name"] == "My category"


def test_bar_plot_cats_dicts_with_typo():
    """
    A typo in the cat properties dict - will fill in the name from the dict key
    """
    plot = _verify_rendered(
        bargraph.plot(
            {"Sample1": {"Cat1": 2}},
            {"Cat1": {"name_with_typo": "My category"}},
            {"id": "test_bar_plot_cats_dicts_with_typo", "title": "Test: Bar Graph"},
        )
    )

    assert report.plot_data[plot.anchor]["datasets"][0]["cats"][0]["name"] == "Cat1"


def test_bar_plot_cats_mismatch_cats_and_ds_count():
    """
    Multiple datasets, but the list lengths are not matching between data and cats. Should throw an error
    """
    with pytest.raises(RunError):
        bargraph.plot(
            [{"Sample1": {"Cat1": 2}}],
            [{"Cat1": {"name": "My category"}}, {"Cat1": {"name": "My category"}}],
            {"id": "test_bar_plot_cats_mismatch_cats_and_ds_count", "title": "Test: Bar Graph"},
        )


def test_bar_plot_fill_cats():
    """
    Multiple datasets, but only one dict of cats - in this case, should copy these cats between datasets
    """
    plot = _verify_rendered(
        bargraph.plot(
            data=[{"Sample1": {"Cat1": 2, "Cat2": 2}}, {"Sample1": {"Cat1": 2, "Cat3": 2}}],
            cats={"Cat1": {"name": "My category"}},
            pconfig={"id": "test_bar_plot_fill_cats", "title": "Test: Bar Graph"},
        )
    )
    assert len(report.plot_data[plot.anchor]["datasets"]) == 2
    assert len(report.plot_data[plot.anchor]["datasets"][0]["cats"]) == 1
    assert len(report.plot_data[plot.anchor]["datasets"][1]["cats"]) == 1
    assert report.plot_data[plot.anchor]["datasets"][0]["cats"][0]["name"] == "My category"
    assert report.plot_data[plot.anchor]["datasets"][1]["cats"][0]["name"] == "My category"


def test_bar_plot_no_cats():
    """
    Cats parameter is missing, will fill in one from the data
    """

    plot = _verify_rendered(
        bargraph.plot(
            {
                "Sample1": {"Cat1": 2, "Cat2": 2},
                "Sample2": {"Cat1": 1, "Cat3": 1},
            },
            pconfig={"id": "test_bar_plot_no_cats", "title": "Test: Bar Graph"},
        )
    )

    assert len(report.plot_data[plot.anchor]["datasets"][0]["cats"]) == 3


def test_bar_plot_sample_groups():
    """
    Test sample_groups configuration for visual grouping
    """
    plot = _verify_rendered(
        bargraph.plot(
            {
                "Sample1": {"Cat1": 10, "Cat2": 20},
                "Sample2": {"Cat1": 15, "Cat2": 25},
                "Sample3": {"Cat1": 12, "Cat2": 22},
                "Sample4": {"Cat1": 18, "Cat2": 28},
            },
            ["Cat1", "Cat2"],
            {
                "id": "test_bar_plot_sample_groups",
                "title": "Test: Bar Graph with Sample Groups",
                "sample_groups": {
                    "Group 1": [["Sample1", "Sample1"], ["Sample2", "Sample2"]],
                    "Group 2": [["Sample3", "Sample3"], ["Sample4", "Sample4"]],
                },
            },
        )
    )

    ds = report.plot_data[plot.anchor]["datasets"][0]
    # Samples should be reordered according to groups
    assert ds["samples"] == ["Sample4", "Sample3", "Sample2", "Sample1"]  # reversed for display
    # Group labels should be present
    assert ds["group_labels"] == ["Group 2", "Group 2", "Group 1", "Group 1"]  # reversed


def test_bar_plot_sample_groups_with_names():
    """
    Test sample_groups with custom group names (now directly in dict keys)
    """
    plot = _verify_rendered(
        bargraph.plot(
            {
                "Sample1": {"Cat1": 10},
                "Sample2": {"Cat1": 15},
            },
            ["Cat1"],
            {
                "id": "test_bar_plot_sample_groups_with_names",
                "title": "Test: Bar Graph with Named Groups",
                "sample_groups": {
                    "Condition A": [["Sample1", "Sample1"]],
                    "Condition B": [["Sample2", "Sample2"]],
                },
            },
        )
    )

    ds = report.plot_data[plot.anchor]["datasets"][0]
    # Custom group names should be used
    assert ds["group_labels"] == ["Condition B", "Condition A"]  # reversed


def test_bar_plot_sample_groups_ungrouped():
    """
    Test that samples not in any group get added to 'Other'
    """
    plot = _verify_rendered(
        bargraph.plot(
            {
                "Sample1": {"Cat1": 10},
                "Sample2": {"Cat1": 15},
                "Sample3": {"Cat1": 20},
            },
            ["Cat1"],
            {
                "id": "test_bar_plot_sample_groups_ungrouped",
                "title": "Test: Bar Graph with Ungrouped Samples",
                "sample_groups": {"Group 1": [["Sample1", "Sample1"]]},  # Sample2 and Sample3 not in any group
            },
        )
    )

    ds = report.plot_data[plot.anchor]["datasets"][0]
    # Sample1 should be in Group 1, others in Other
    # Order: grouped samples first (Sample1), then ungrouped (Sample2, Sample3)
    assert "Other" in ds["group_labels"]
    assert "Group 1" in ds["group_labels"]


def test_bar_plot_sample_groups_disables_sort():
    """
    Test that sample_groups disables sort_samples
    """
    inputs = bargraph.BarPlotInputData.create(
        {"Sample1": {"Cat1": 10}, "Sample2": {"Cat1": 15}},
        ["Cat1"],
        {
            "id": "test_bar_plot_sample_groups_disables_sort",
            "title": "Test",
            "sample_groups": {"Group 1": [["Sample1", "Sample1"]], "Group 2": [["Sample2", "Sample2"]]},
            "sort_samples": True,  # Should be overridden
        },
    )

    assert inputs.pconfig.sort_samples is False


def test_bar_plot_sample_groups_disables_clustering():
    """
    Test that sample_groups disables cluster_samples
    """
    inputs = bargraph.BarPlotInputData.create(
        {"Sample1": {"Cat1": 10}, "Sample2": {"Cat1": 15}},
        ["Cat1"],
        {
            "id": "test_bar_plot_sample_groups_disables_clustering",
            "title": "Test",
            "sample_groups": {"Group 1": [["Sample1", "Sample1"]], "Group 2": [["Sample2", "Sample2"]]},
            "cluster_samples": True,  # Should be overridden
        },
    )

    assert inputs.pconfig.cluster_samples is False


def test_bar_plot_sample_groups_empty_group():
    """
    Test that empty groups (groups with no matching samples) are handled gracefully
    """
    plot = _verify_rendered(
        bargraph.plot(
            {
                "Sample1": {"Cat1": 10},
                "Sample2": {"Cat1": 15},
            },
            ["Cat1"],
            {
                "id": "test_bar_plot_sample_groups_empty_group",
                "title": "Test: Bar Graph with Empty Group",
                "sample_groups": {
                    "Group A": [["Sample1", "Sample1"]],
                    "Empty Group": [["NonExistentSample", "NonExistent"]],  # This group has no matching samples
                    "Group B": [["Sample2", "Sample2"]],
                },
            },
        )
    )

    ds = report.plot_data[plot.anchor]["datasets"][0]
    # Only samples that exist should be in the output
    # Empty group should not contribute any samples or labels
    assert len(ds["samples"]) == 2
    assert len(ds["group_labels"]) == 2
    # Group labels should be "Group A" and "Group B" (no "Empty Group")
    assert "Empty Group" not in ds["group_labels"]


def test_bar_plot_sample_groups_multiple_entries():
    """
    Test same sample appearing in multiple groups with lists for offset alignment
    """
    plot = _verify_rendered(
        bargraph.plot(
            {
                "Sample1_25nt": {"Frame0": 50, "Frame1": 30, "Frame2": 20},
                "Sample1_26nt": {"Frame0": 60, "Frame1": 25, "Frame2": 15},
                "Sample2_25nt": {"Frame0": 55, "Frame1": 28, "Frame2": 17},
                "Sample2_26nt": {"Frame0": 65, "Frame1": 22, "Frame2": 13},
            },
            ["Frame0", "Frame1", "Frame2"],
            {
                "id": "test_bar_plot_sample_groups_multiple_entries",
                "title": "Test: Bar Graph with Multiple Entries Per Sample",
                "sample_groups": {
                    "25nt": [["Sample1_25nt", "Sample1"], ["Sample2_25nt", "Sample2"]],
                    "26nt": [["Sample1_26nt", "Sample1"], ["Sample2_26nt", "Sample2"]],
                },
            },
        )
    )

    ds = report.plot_data[plot.anchor]["datasets"][0]
    # All 4 samples should be present
    assert len(ds["samples"]) == 4
    # Group labels should have 2 of each type
    assert ds["group_labels"].count("25nt") == 2
    assert ds["group_labels"].count("26nt") == 2
    # Offset groups should map sample keys to their base sample names
    assert ds["offset_groups"]["Sample1_25nt"] == "Sample1"
    assert ds["offset_groups"]["Sample1_26nt"] == "Sample1"
    assert ds["offset_groups"]["Sample2_25nt"] == "Sample2"
    assert ds["offset_groups"]["Sample2_26nt"] == "Sample2"


def test_linegraph_axis_controlled_by_switches_valid():
    """Test that valid axis_controlled_by_switches values are accepted."""
    # Test with yaxis only (default behavior)
    config1 = LinePlotConfig(id="test1", title="Test", axis_controlled_by_switches=["yaxis"])
    assert config1.axis_controlled_by_switches == ["yaxis"]

    # Test with xaxis only
    config2 = LinePlotConfig(id="test2", title="Test", axis_controlled_by_switches=["xaxis"])
    assert config2.axis_controlled_by_switches == ["xaxis"]

    # Test with both axes
    config3 = LinePlotConfig(id="test3", title="Test", axis_controlled_by_switches=["xaxis", "yaxis"])
    assert config3.axis_controlled_by_switches == ["xaxis", "yaxis"]

    # Test with None (default)
    config4 = LinePlotConfig(id="test4", title="Test")
    assert config4.axis_controlled_by_switches is None


def test_linegraph_axis_controlled_by_switches_invalid():
    """Test that invalid axis_controlled_by_switches values are rejected with a useful error."""
    with patch("logging.Logger.error") as err:
        config = LinePlotConfig(id="test", title="Test", axis_controlled_by_switches=["invalid"])
        assert config.axis_controlled_by_switches is None
        errs = "\n".join(call.args[0] for call in err.mock_calls if call.args)
        assert "'axis_controlled_by_switches'" in errs
        assert "Literal['xaxis', 'yaxis']" in errs


def test_linegraph_axis_controlled_by_switches_string_instead_of_list():
    """Test that a flat string instead of a list is rejected with a useful error."""
    with patch("logging.Logger.error") as err:
        config = LinePlotConfig(id="test", title="Test", axis_controlled_by_switches="yaxis")  # type: ignore
        assert config.axis_controlled_by_switches is None
        errs = "\n".join(call.args[0] for call in err.mock_calls if call.args)
        assert "'axis_controlled_by_switches'" in errs
        assert "List" in errs


def test_linegraph_axis_controlled_by_switches_in_plot():
    """Test that axis_controlled_by_switches works in actual plot creation."""
    dataset = {"Sample1": {0: 1, 1: 2}}

    # Test with xaxis
    plot = _verify_rendered(
        linegraph.plot(
            dataset,
            LinePlotConfig(id="test_axis_xaxis", title="Test", axis_controlled_by_switches=["xaxis"]),
        )
    )
    assert isinstance(plot, linegraph.LinePlot)


def test_linegraph_smooth():
    SMOOTH_TO = 2
    dataset = {"Smoothed": {0: 1, 1: 1, 2: 1}, "Unsmoothed": {0: 1, 1: 1}}

    plot = _verify_rendered(
        linegraph.plot(
            dataset,
            {"id": "test_linegraph_smooth", "title": "Test: Line Graph", "smooth_points": SMOOTH_TO},
        )
    )

    for in_series, out_series in zip(dataset.values(), report.plot_data[plot.anchor]["datasets"][0]["lines"]):
        assert min(len(in_series), SMOOTH_TO) == len(out_series["pairs"])


def test_linegraph_multiple_datasets():
    plot = _verify_rendered(
        linegraph.plot(
            [{"Sample1": {0: 1, 1: 1}}, {"Sample1": {0: 2, 1: 2}}],
            {
                "id": "test_linegraph_multiple_datasets",
                "title": "Test: Line Graph",
                "data_labels": ["Dataset1", "Dataset2"],
            },
        )
    )

    assert len(report.plot_data[plot.anchor]["datasets"]) == 2


@pytest.mark.parametrize(
    "development,export_plots,export_plot_formats",
    [
        (False, False, None),  # default mode - embed, no export
        (False, True, None),  # embed + export all formats
        (True, True, ["pdf"]),  # link png + export pdf (should also export png for html)
        (True, False, None),  # link png + no export (should only export png)
    ],
)
@pytest.mark.filterwarnings("ignore:setDaemon")
@pytest.mark.skip(reason="Fails on CI")
def test_flat_plot(tmp_path, monkeypatch, development, export_plot_formats, export_plots):
    monkeypatch.setattr(tempfile, "mkdtemp", lambda *args, **kwargs: tmp_path)

    plot_id = "test_plot"
    plot = linegraph.plot(
        {"Sample1": {0: 1, 1: 1}},
        {"id": plot_id, "title": "Line Graph"},
    )
    assert isinstance(plot, Plot)

    plot.flat = True
    config.development = development
    config.export_plots = export_plots
    if export_plot_formats:
        config.export_plot_formats = export_plot_formats

    html = plot.add_to_report(
        module_anchor=Anchor("test"), section_anchor=Anchor("test"), plots_dir_name=config.plots_dir_name
    )
    # Process any batched exports
    process_batch_exports()

    assert len(report.plot_data) == 0
    assert html is not None
    if not development:
        assert f'<div class="mqc_mplplot" style="" id="{plot_id}"><img src="data:image/png;base64' in html
        if not export_plots:
            for fmt in ["png", "pdf", "svg"]:
                assert not (tmp_path / f"multiqc_plots/{fmt}/{plot_id}.{fmt}").is_file()
    else:
        assert f'<div class="mqc_mplplot" style="" id="{plot_id}"><img src="multiqc_plots/png/{plot_id}.png' in html
        assert (tmp_path / f"multiqc_plots/png/{plot_id}.png").is_file()
        assert (tmp_path / f"multiqc_plots/png/{plot_id}.png").stat().st_size > 0
        if not export_plots:
            for fmt in ["pdf", "svg"]:
                assert not (tmp_path / f"multiqc_plots/{fmt}/{plot_id}.{fmt}").is_file()
    if export_plots:
        for fmt in export_plot_formats or ["png", "pdf", "svg"]:
            assert (tmp_path / f"multiqc_plots/{fmt}/{plot_id}.{fmt}").is_file()
            assert (tmp_path / f"multiqc_plots/{fmt}/{plot_id}.{fmt}").stat().st_size > 0


def test_missing_pconfig(reset):
    from multiqc import config

    config.strict = True

    linegraph.plot({"Sample1": {0: 1, 1: 1}})
    assert report.lint_errors == [
        "pconfig with required fields 'id' and 'title' must be provided for plot LinePlotConfig",
    ]

    _verify_rendered(linegraph.plot({"Sample1": {0: 1, 1: 1}}))
    plot_id = list(report.plot_data.keys())[0]
    assert plot_id.startswith("lineplot-")


@pytest.mark.parametrize("strict", [True, False])
def test_incorrect_fields(strict, reset):
    from multiqc import config

    config.strict = strict

    pconfig = {
        "id": "test_incorrect_fields",
        "title": "Test: Line Graph",
        "unknown_field": "value",
        "x_lines": "wrong_type",
    }

    if strict:
        with pytest.raises(ModuleConfigValidationError):
            linegraph.plot({"Sample1": {0: 1, 1: 1}}, pconfig=pconfig)
    else:
        with patch("logging.Logger.error") as err, patch("logging.Logger.warning") as warn:
            _verify_rendered(linegraph.plot({"Sample1": {0: 1, 1: 1}}, pconfig=pconfig))
            errs = "\n".join(call.args[0] for call in err.mock_calls if call.args)
            assert "• 'x_lines': failed to parse value 'wrong_type'" in errs
            assert "errors while parsing lineplot.pconfig[id='test_incorrect_fields']" in errs
            warnings = "\n".join(call.args[0] for call in warn.mock_calls if call.args)
            assert "• 'unknown_field': unrecognized field" in warnings
        assert "test_incorrect_fields" in report.plot_data


@pytest.mark.parametrize("strict", [True, False])
def test_missing_id_and_title(strict, reset):
    from multiqc import config

    config.strict = strict
    if strict:
        with pytest.raises(ModuleConfigValidationError):
            linegraph.plot({"Sample1": {0: 1, 1: 1}}, pconfig={})
    else:
        with patch("logging.Logger.error") as log:
            _verify_rendered(linegraph.plot({"Sample1": {0: 1, 1: 1}}, pconfig={}))
            errs = "\n".join(call.args[0] for call in log.mock_calls if call.args)
            assert "• 'id': missing required field" in errs
            assert "• 'title': missing required field" in errs
        plot_id = list(report.plot_data.keys())[0]
        assert plot_id.startswith("lineplot-")


def test_incorrect_color():
    with patch("logging.Logger.error") as err:
        _verify_rendered(
            linegraph.plot(
                {"Sample1": {0: 1, 1: 1}},
                pconfig={
                    "id": "test_incorrect_color",
                    "title": "Line Graph",
                    "extra_series": [{"color": "invalid"}],
                },
            )
        )
        errs = "\n".join(call.args[0] for call in err.mock_calls if call.args)
        assert "• 'color': invalid color value 'invalid'" in errs


def test_extra_series_multiple_datasets():
    """Should zip series with the datasets"""
    plot_id = "my_plot"
    _verify_rendered(
        linegraph.plot(
            [{"Sample1": {0: 1, 1: 1}}, {"Sample1": {0: 2, 1: 2}}],
            pconfig=LinePlotConfig(
                id=plot_id,
                title="Line Graph",
                extra_series=[Series(pairs=[(1, 2)], name="Extra1")],
            ),
        )
    )

    anchor = Anchor(plot_id)
    assert len(report.plot_data[anchor]["datasets"][0]["lines"]) == 2
    assert len(report.plot_data[anchor]["datasets"][0]["lines"][0]["pairs"]) == 2
    assert len(report.plot_data[anchor]["datasets"][0]["lines"][1]["pairs"]) == 1
    assert report.plot_data[anchor]["datasets"][0]["lines"][0]["name"] == "Sample1"
    assert report.plot_data[anchor]["datasets"][0]["lines"][1]["name"] == "Extra1"

    assert len(report.plot_data[anchor]["datasets"][1]["lines"]) == 2
    assert len(report.plot_data[anchor]["datasets"][1]["lines"][0]["pairs"]) == 2
    assert len(report.plot_data[anchor]["datasets"][1]["lines"][1]["pairs"]) == 1
    assert report.plot_data[anchor]["datasets"][1]["lines"][0]["name"] == "Sample1"
    assert report.plot_data[anchor]["datasets"][1]["lines"][1]["name"] == "Extra1"


def test_multiple_extra_series():
    """Should add two series to the dataset in addition to the main data"""
    plot_id = "my_plot"
    _verify_rendered(
        linegraph.plot(
            {"Sample1": {0: 1, 1: 1}},
            pconfig=LinePlotConfig(
                id=plot_id,
                title="Line Graph",
                extra_series=[Series(pairs=[(1, 2)], name="Extra1"), Series(pairs=[(2, 3)], name="Extra2")],
            ),
        )
    )

    anchor = Anchor(plot_id)
    assert len(report.plot_data[anchor]["datasets"]) == 1
    assert len(report.plot_data[anchor]["datasets"][0]["lines"]) == 3
    assert len(report.plot_data[anchor]["datasets"][0]["lines"][0]["pairs"]) == 2
    assert len(report.plot_data[anchor]["datasets"][0]["lines"][1]["pairs"]) == 1
    assert len(report.plot_data[anchor]["datasets"][0]["lines"][2]["pairs"]) == 1
    assert report.plot_data[anchor]["datasets"][0]["lines"][0]["name"] == "Sample1"
    assert report.plot_data[anchor]["datasets"][0]["lines"][1]["name"] == "Extra1"
    assert report.plot_data[anchor]["datasets"][0]["lines"][2]["name"] == "Extra2"


def test_extra_series_multiple_datasets_different_series():
    """Should zip series with the datasets"""
    plot_id = "my_plot"
    datasets = [{"Sample1": {0: 1, 1: 1}}, {"Sample1": {0: 2, 1: 2}}]
    _verify_rendered(
        linegraph.plot(
            datasets,
            pconfig=LinePlotConfig(
                id=plot_id,
                title="Line Graph",
                extra_series=[[Series(pairs=[(1, 2)], name="Extra1")], [Series(pairs=[(2, 3)], name="Extra2")]],
            ),
        )
    )

    anchor = Anchor(plot_id)
    assert len(report.plot_data[anchor]["datasets"]) == 2
    for ds in report.plot_data[anchor]["datasets"]:
        assert len(ds["lines"]) == 2
        assert len(ds["lines"][0]["pairs"]) == 2
        assert len(ds["lines"][1]["pairs"]) == 1
        assert ds["lines"][0]["name"] == "Sample1"
    assert report.plot_data[anchor]["datasets"][0]["lines"][1]["name"] == "Extra1"
    assert report.plot_data[anchor]["datasets"][1]["lines"][1]["name"] == "Extra2"


def test_extra_series_multiple_datasets_multiple_series():
    """Should copy the extra series to all datasets"""
    plot_id = "my_plot"
    _verify_rendered(
        linegraph.plot(
            [{"Sample1": {0: 1, 1: 1}}, {"Sample1": {0: 2, 1: 2}}],
            pconfig=LinePlotConfig(
                id=plot_id,
                title="Line Graph",
                extra_series=[Series(pairs=[(1, 2)], name="Extra1"), Series(pairs=[(2, 3)], name="Extra2")],
            ),
        )
    )

    anchor = Anchor(plot_id)
    assert len(report.plot_data[anchor]["datasets"]) == 2
    for ds in report.plot_data[anchor]["datasets"]:
        assert len(ds["lines"]) == 3
        assert len(ds["lines"][0]["pairs"]) == 2
        assert len(ds["lines"][1]["pairs"]) == 1
        assert len(ds["lines"][2]["pairs"]) == 1
        assert ds["lines"][0]["name"] == "Sample1"
        assert ds["lines"][1]["name"] == "Extra1"
        assert ds["lines"][2]["name"] == "Extra2"


def test_dash_styles():
    plot_id = "my_plot"
    pconfig = {
        "id": plot_id,
        "title": "Line Graph",
        "extra_series": [
            {"dash": "dash", "pairs": [(1, 1)]},
            {"dashStyle": "dash", "pairs": [(1, 1)]},
            {"dash": "ShortDash", "pairs": [(1, 1)]},
            {"dashStyle": "ShortDash", "pairs": [(1, 1)]},
        ],
    }
    data = {
        "Sample1": {0: 1, 1: 1},
    }
    anchor = Anchor(plot_id)
    with patch("logging.Logger.warning") as log:
        _verify_rendered(linegraph.plot(data, pconfig=pconfig))
        warnings = "\n".join(call.args[0] for call in log.mock_calls if call.args)
        assert "• 'extra_series': 'dashStyle' field is deprecated. Please use 'dash' instead" in warnings
        assert "• 'dash': 'ShortDash' is a deprecated dash style, use 'dash'" in warnings
    assert len(report.plot_data[anchor]["datasets"][0]["lines"]) == 5
    for line in report.plot_data[anchor]["datasets"][0]["lines"][1:]:
        assert line["dash"] == "dash"


def test_table_default_sort():
    from multiqc.plots.table_object import _get_sortlist_js

    headers: Dict[str, ColumnDict] = {"x": {"title": "Metric X"}, "y": {"title": "Metric Y"}}
    p = table.plot(
        data={
            "sample1": {"x": 1, "y": 2},
            "sample2": {"x": 3, "y": 4},
        },
        headers=headers,
        pconfig=table.TableConfig(
            id="table",
            title="Table",
            defaultsort=[
                {"column": "y", "direction": "desc"},
                {"column": "x", "direction": "asc"},
            ],
        ),
    )
    assert isinstance(p, Plot)
    sort_string = _get_sortlist_js(p.datasets[0].dt)
    assert sort_string == "[[2, 1], [1, 0]]"


def test_table_custom_plot_config_hidden(reset):
    """
    Test that custom_plot_config can set column properties at the table level.
    When 'hidden: true' is set at the table level, all columns should be hidden.
    """
    table_id = "test_table_hidden"

    # Set custom_plot_config for this table
    config.custom_plot_config = {
        table_id: {
            "hidden": True,  # Should apply to all columns
        }
    }

    headers: Dict[str, ColumnDict] = {
        "x": {"title": "Metric X"},
        "y": {"title": "Metric Y"},
        "z": {"title": "Metric Z"},
    }

    p = table.plot(
        data={
            "sample1": {"x": 1, "y": 2, "z": 3},
            "sample2": {"x": 4, "y": 5, "z": 6},
        },
        headers=headers,
        pconfig=table.TableConfig(id=table_id, title="Test Table"),
    )

    assert isinstance(p, Plot)

    # Check that all columns are hidden
    dt = p.datasets[0].dt
    for section in dt.section_by_id.values():
        for col_key, col_meta in section.column_by_key.items():
            assert col_meta.hidden is True, f"Column {col_key} should be hidden"


def test_table_custom_plot_config_scale(reset):
    """
    Test that custom_plot_config can set the color scale at the table level.
    When 'scale: RdYlGn' is set at the table level, all columns should use that scale.
    """
    table_id = "test_table_scale"

    # Set custom_plot_config for this table
    config.custom_plot_config = {
        table_id: {
            "scale": "RdYlGn",  # Should apply to all columns
        }
    }

    headers: Dict[str, ColumnDict] = {
        "x": {"title": "Metric X", "scale": "Blues"},  # This should be overridden
        "y": {"title": "Metric Y", "scale": "Reds"},  # This should be overridden
        "z": {"title": "Metric Z"},  # This should get RdYlGn
    }

    p = table.plot(
        data={
            "sample1": {"x": 1, "y": 2, "z": 3},
            "sample2": {"x": 4, "y": 5, "z": 6},
        },
        headers=headers,
        pconfig=table.TableConfig(id=table_id, title="Test Table"),
    )

    assert isinstance(p, Plot)

    # Check that all columns have the RdYlGn scale
    dt = p.datasets[0].dt
    for section in dt.section_by_id.values():
        for col_key, col_meta in section.column_by_key.items():
            assert col_meta.scale == "RdYlGn", f"Column {col_key} should have scale 'RdYlGn', got '{col_meta.scale}'"


def test_table_custom_plot_config_multiple_properties(reset):
    """
    Test that custom_plot_config can set multiple column properties at once.
    """
    table_id = "test_table_multi"

    # Set multiple properties at the table level
    config.custom_plot_config = {
        table_id: {
            "hidden": False,
            "scale": "Purples",
            "suffix": " units",
        }
    }

    headers: Dict[str, ColumnDict] = {
        "x": {"title": "Metric X", "hidden": True},  # Should be overridden to False
        "y": {"title": "Metric Y"},
    }

    p = table.plot(
        data={
            "sample1": {"x": 1, "y": 2},
            "sample2": {"x": 3, "y": 4},
        },
        headers=headers,
        pconfig=table.TableConfig(id=table_id, title="Test Table"),
    )

    assert isinstance(p, Plot)

    # Check that all properties are applied
    dt = p.datasets[0].dt
    for section in dt.section_by_id.values():
        for col_key, col_meta in section.column_by_key.items():
            assert col_meta.hidden is False, f"Column {col_key} should not be hidden"
            assert col_meta.scale == "Purples", f"Column {col_key} should have scale 'Purples'"
            assert col_meta.suffix == " units", f"Column {col_key} should have suffix ' units'"


def test_table_custom_plot_config_invalid_field(reset):
    """
    Test that invalid fields in custom_plot_config are silently ignored.
    This should not crash when a table-level property doesn't exist on TableConfig.
    """
    table_id = "test_table_invalid"

    # Set invalid properties - 'hidden' is not a TableConfig field, only a ColumnMeta field
    config.custom_plot_config = {
        table_id: {
            "hidden": True,  # Valid ColumnMeta field, should apply to columns
            "invalid_field": "value",  # Invalid field, should be ignored
        }
    }

    headers: Dict[str, ColumnDict] = {
        "x": {"title": "Metric X"},
        "y": {"title": "Metric Y"},
    }

    # This should not raise an error
    p = table.plot(
        data={
            "sample1": {"x": 1, "y": 2},
            "sample2": {"x": 3, "y": 4},
        },
        headers=headers,
        pconfig=table.TableConfig(id=table_id, title="Test Table"),
    )

    assert isinstance(p, Plot)

    # Check that the valid field (hidden) was applied
    dt = p.datasets[0].dt
    for section in dt.section_by_id.values():
        for col_key, col_meta in section.column_by_key.items():
            assert col_meta.hidden is True, f"Column {col_key} should be hidden"
