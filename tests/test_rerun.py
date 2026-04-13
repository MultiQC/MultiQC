import difflib
import json
from datetime import datetime, timedelta

import polars as pl

import multiqc
from multiqc.core.update_config import ClConfig
from multiqc.plots.bargraph import BarPlotConfig, BarPlotInputData, CatConf
from multiqc.plots.linegraph import LinePlotConfig, LinePlotNormalizedInputData, Series
from multiqc.plots.plot import PlotType, plot_anchor
from multiqc.types import SampleName


def test_rerun_parquet(data_dir, tmp_path):
    """Test running MultiQC on intermediate output multiqc_data.parquet.
    Run MultiQC on a set of inputs.
    Then run MultiQC on the multiqc_data.parquet from the first run. The reports should be identical.
    It should only work with multiqc_data.parquet is passed explicitly. Otherwise MultiQC output folder should be ignored.
    """
    # Run 1: Normal run on test data
    run_a_dir = tmp_path / "run_a"
    run_a_dir.mkdir()
    multiqc.run(data_dir / "modules/fastp/SAMPLE.json", cfg=ClConfig(output_dir=run_a_dir, strict=True))

    # Get the first report contents
    with open(run_a_dir / "multiqc_data" / "multiqc_data.json") as f:
        report1_data = json.load(f)

    # Run 2: Run on the intermediate data from run1
    run_b_dir = tmp_path / "run_b"
    run_b_dir.mkdir()
    multiqc.run(run_a_dir / "multiqc_data" / "multiqc.parquet", cfg=ClConfig(output_dir=run_b_dir, strict=True))

    # Compare reports
    with open(run_b_dir / "multiqc_data" / "multiqc_data.json") as f:
        report2_data = json.load(f)

    # Compare only the relevant fields, not the entire report
    # The 'config_analysis_dir_abs' will be different between runs
    for key in ["report_plot_data"]:
        assert key in report1_data, f"Key {key} missing from first report"
        assert key in report2_data, f"Key {key} missing from second report"
        assert report1_data[key] == report2_data[key], f"Value for {key} differs between reports"


def test_rerun_and_combine(data_dir, tmp_path):
    """Test adding new data to report.
    Run MultiQC on a set of inputs.
    Then run MultiQC on the intermediate output multiqc_data.json of that run, plus add new inputs.
    The output should be the same as if MultiQC was run on the combined inputs directly.
    """
    # Run MultiQC on inputs A
    run_a_dir = tmp_path / "run_a"
    multiqc.run(
        data_dir / "modules/fastp/single_end",
        cfg=ClConfig(
            output_dir=run_a_dir,
            strict=True,
        ),
    )

    # Run MultiQC on outputs of A + inputs B
    run_combined_dir = tmp_path / "run_combined"
    multiqc.run(
        data_dir / "modules/fastp/SAMPLE.json",
        run_a_dir / "multiqc_data" / "multiqc.parquet",
        cfg=ClConfig(output_dir=run_combined_dir, strict=True),
    )

    with open(run_combined_dir / "multiqc_data" / "multiqc_data.json") as f:
        report_combined_data = json.load(f)

    # Run MultiQC on inputs A + B directly
    run_normal_dir = tmp_path / "run_normal"
    multiqc.run(
        data_dir / "modules/fastp/SAMPLE.json",  # "smalltest_S10_R1_001"
        data_dir / "modules/fastp/single_end",  # "SRR5442949_1"
        cfg=ClConfig(output_dir=run_normal_dir, strict=True),
    )

    with open(run_normal_dir / "multiqc_data" / "multiqc_data.json") as f:
        report_normal_data = json.load(f)

    # Write report plot data to separate files for inspection
    combined_plot_data_file = tmp_path / "combined_plot_data.json"
    normal_plot_data_file = tmp_path / "normal_plot_data.json"
    with open(combined_plot_data_file, "w") as f:
        json.dump(report_combined_data["report_plot_data"], f, indent=4)
    with open(normal_plot_data_file, "w") as f:
        json.dump(report_normal_data["report_plot_data"], f, indent=4)
    print(f"\nCombined plot data written to: {combined_plot_data_file}")
    print(f"Normal plot data written to: {normal_plot_data_file}")

    # Compare only the relevant fields, not the entire report
    # The 'config_analysis_dir_abs' will be different between runs
    for key in [
        "report_data_sources",
        "report_plot_data",
    ]:
        assert key in report_combined_data, f"Key {key} missing from combined report"
        assert key in report_normal_data, f"Key {key} missing from direct report"
        assert report_combined_data[key] == report_normal_data[key], f"Value for {key} differs between reports"


def test_merge_linegraph():
    """Test merging two linegraph inputs.
    Create two different datasets with some overlapping samples.
    Merge them and verify that the overlapping sample from the first dataset is replaced
    by the one from the second dataset, and non-overlapping samples are preserved.
    """
    # Create plot config
    plot_id = "test_merge_plot"
    pconfig = LinePlotConfig(id=plot_id, title="Test Merge Plot")
    anchor = plot_anchor(pconfig)

    # Create first input data - two samples: "Sample1" and "Sample2"
    dataset1 = [
        Series(name="Sample1", pairs=[(0, 1), (1, 2), (2, 3)]),
        Series(name="Sample2", pairs=[(0, 2), (1, 3), (2, 4)]),
    ]
    input_data1 = LinePlotNormalizedInputData(
        anchor=anchor,
        plot_type=PlotType.LINE,
        data=[dataset1],
        pconfig=pconfig,
        sample_names=[SampleName("Sample1"), SampleName("Sample2")],
        creation_date=datetime.now() - timedelta(days=1),
    )

    # Create second input data - two samples: "Sample1" (overlapping) and "Sample3" (new)
    dataset2 = [
        Series(name="Sample1", pairs=[(0, 5), (1, 6), (2, 7)]),  # Different data for Sample1
        Series(name="Sample3", pairs=[(0, 3), (1, 4), (2, 5)]),  # New sample
    ]
    input_data2 = LinePlotNormalizedInputData(
        anchor=anchor,
        plot_type=PlotType.LINE,
        data=[dataset2],
        pconfig=pconfig,
        sample_names=[SampleName("Sample1"), SampleName("Sample3")],
        creation_date=datetime.now(),
    )

    # Merge the two inputs
    merged_data = LinePlotNormalizedInputData.merge(input_data1, input_data2)

    # Convert to dataframe for verification - this helps validate the core merging logic
    merged_df = merged_data.to_df()

    # Verify the merged data has three unique samples: Sample1 (from dataset2), Sample2, and Sample3
    unique_samples = merged_df.select("sample").unique().to_series()
    assert len(unique_samples) == 3
    assert "Sample1" in unique_samples
    assert "Sample2" in unique_samples
    assert "Sample3" in unique_samples

    # Group by sample and verify data points
    sample1_data = merged_df.filter(pl.col("sample") == "Sample1")
    sample2_data = merged_df.filter(pl.col("sample") == "Sample2")
    sample3_data = merged_df.filter(pl.col("sample") == "Sample3")

    # Sample1 should have the values from dataset2
    assert sample1_data.height == 3  # 3 data points
    sample1_y_vals = sample1_data.select("y_val").sort("y_val").to_series().to_list()
    assert sample1_y_vals == ["5", "6", "7"]  # Values come from dataset2

    # Sample2 should have values from dataset1 (unchanged)
    assert sample2_data.height == 3
    sample2_y_vals = sample2_data.select("y_val").sort("y_val").to_series().to_list()
    assert sample2_y_vals == ["2", "3", "4"]

    # Sample3 should have values from dataset2
    assert sample3_data.height == 3
    sample3_y_vals = sample3_data.select("y_val").sort("y_val").to_series().to_list()
    assert sample3_y_vals == ["3", "4", "5"]


def test_merge_bargraph():
    """Test merging two bar graph inputs.
    Create two different datasets with some overlapping samples and categories.
    Merge them and verify that the overlapping sample/category data from the first dataset
    is replaced by data from the second dataset, and non-overlapping samples are preserved.
    """
    # Create plot config
    plot_id = "test_bargraph_merge"
    pconfig = BarPlotConfig(id=plot_id, title="Test Bar Graph Merge")
    anchor = plot_anchor(pconfig)

    # Create first input data - two samples with two categories each
    dataset1 = {
        "Sample1": {"Cat1": 10, "Cat2": 15},
        "Sample2": {"Cat1": 20, "Cat2": 25},
    }

    cats1 = {
        "Cat1": CatConf(name="Category 1", color="#ff0000"),
        "Cat2": CatConf(name="Category 2", color="#00ff00"),
    }

    input_data1 = BarPlotInputData(
        anchor=anchor,
        plot_type=PlotType.BAR,
        data=[dataset1],  # type: ignore
        cats=[cats1],  # type: ignore
        pconfig=pconfig,
        creation_date=datetime.now() - timedelta(days=1),
    )

    # Create second input data - two samples: "Sample1" (overlapping) and "Sample3" (new)
    # with two categories: "Cat1" (overlapping) and "Cat3" (new)
    dataset2 = {
        "Sample1": {"Cat1": 30, "Cat3": 35},  # Different data for Sample1/Cat1 and new category
        "Sample3": {"Cat1": 40, "Cat3": 45},  # New sample
    }

    cats2 = {
        "Cat1": CatConf(name="Category 1", color="#ff0000"),
        "Cat3": CatConf(name="Category 3", color="#0000ff"),
    }

    input_data2 = BarPlotInputData(
        anchor=anchor,
        plot_type=PlotType.BAR,
        data=[dataset2],  # type: ignore
        cats=[cats2],  # type: ignore
        pconfig=pconfig,
        creation_date=datetime.now(),
    )

    # Merge the two inputs
    merged_data = BarPlotInputData.merge(input_data1, input_data2)

    # Convert to dataframe for verification
    merged_df = merged_data.to_df()

    # Verify the merged data has all samples and categories preserved correctly
    unique_samples = merged_df.select("sample").unique().to_series()
    unique_cats = merged_df.select("category").unique().to_series()

    assert len(unique_samples) == 3
    assert "Sample1" in unique_samples
    assert "Sample2" in unique_samples
    assert "Sample3" in unique_samples

    assert len(unique_cats) == 3
    assert "Cat1" in unique_cats
    assert "Cat2" in unique_cats
    assert "Cat3" in unique_cats

    # Group by sample/category and verify data points
    # Sample1/Cat1 should have the value from dataset2
    sample1_cat1 = merged_df.filter((pl.col("sample") == "Sample1") & (pl.col("category") == "Cat1"))
    assert sample1_cat1.height == 1
    assert float(sample1_cat1.select("bar_value").item()) == 30.0  # Updated value from dataset2

    # Sample1/Cat2 should be overridden
    sample1_cat2 = merged_df.filter((pl.col("sample") == "Sample1") & (pl.col("category") == "Cat2"))
    assert sample1_cat2.height == 0

    # Sample1/Cat3 should be from dataset2
    sample1_cat3 = merged_df.filter((pl.col("sample") == "Sample1") & (pl.col("category") == "Cat3"))
    assert sample1_cat3.height == 1
    assert float(sample1_cat3.select("bar_value").item()) == 35.0

    # Sample3/Cat1 should be from dataset2
    sample3_cat1 = merged_df.filter((pl.col("sample") == "Sample3") & (pl.col("category") == "Cat1"))
    assert sample3_cat1.height == 1
    assert float(sample3_cat1.select("bar_value").item()) == 40.0

    # Sample3/Cat3 should be from dataset2
    sample3_cat3 = merged_df.filter((pl.col("sample") == "Sample3") & (pl.col("category") == "Cat3"))
    assert sample3_cat3.height == 1
    assert float(sample3_cat3.select("bar_value").item()) == 45.0
