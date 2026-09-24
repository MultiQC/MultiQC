import difflib
import json
from datetime import datetime, timedelta

import polars as pl

import multiqc
from multiqc.core import tmp_dir
from multiqc.core.plot_data_store import flush_to_parquet
from multiqc.core.special_case_modules.load_multiqc_data import create_plot_input_data_only
from multiqc.core.update_config import ClConfig
from multiqc.plots import table_object
from multiqc.plots.bargraph import BarPlotConfig, BarPlotInputData, CatConf
from multiqc.plots.linegraph import LinePlotConfig, LinePlotNormalizedInputData, Series
from multiqc.plots.plot import PlotType, plot_anchor
from multiqc.plots.table_object import InputRow, TableConfig
from multiqc.plots.violin import ViolinPlotInputData
from multiqc.types import Anchor, ColumnKey, SampleGroup, SampleName, SectionKey


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


def _make_grouped_violin_data() -> ViolinPlotInputData:
    """
    Build a ViolinPlotInputData with a multi-row group (aggregate + members),
    simulating what table_sample_merge produces.
    """
    pconfig = TableConfig(id="test_grouped", title="Test Grouped Table")
    anchor = plot_anchor(pconfig)
    table_anchor = Anchor(f"{anchor}_table")

    grouped_data: dict = {
        SectionKey("section1"): {
            SampleGroup("SampleA"): [
                InputRow(sample=SampleName("SampleA"), data={"reads": 20000, "quality": 35.0}),
                InputRow(sample=SampleName("SampleA R1"), data={"reads": 10000, "quality": 34.0}),
                InputRow(sample=SampleName("SampleA R2"), data={"reads": 10000, "quality": 36.0}),
            ],
            SampleGroup("SampleB"): [
                InputRow(sample=SampleName("SampleB"), data={"reads": 15000}),
            ],
        }
    }
    headers: dict = {
        SectionKey("section1"): {
            ColumnKey("reads"): {"title": "Total Reads"},
            ColumnKey("quality"): {"title": "Quality"},
        }
    }

    dt = table_object.DataTable.create(
        data=grouped_data,
        table_id=pconfig.id,
        table_anchor=table_anchor,
        pconfig=pconfig.model_copy(),
        headers=headers,
    )

    return ViolinPlotInputData(
        dt=dt,
        plot_type=PlotType.VIOLIN,
        pconfig=pconfig,
        anchor=anchor,
        show_table_by_default=True,
        creation_date=datetime.now(),
    )


def test_grouped_samples_survive_df_roundtrip():
    """
    Verify that sample groups (aggregate row + member rows) survive
    serialization to DataFrame and deserialization back.
    """
    original = _make_grouped_violin_data()

    section = list(original.dt.section_by_id.values())[0]
    assert SampleGroup("SampleA") in section.rows_by_sgroup
    assert len(section.rows_by_sgroup[SampleGroup("SampleA")]) == 3
    assert section.rows_by_sgroup[SampleGroup("SampleA")][0].sample == SampleName("SampleA")
    assert section.rows_by_sgroup[SampleGroup("SampleA")][0].data[ColumnKey("reads")].raw == 20000

    df = original.to_df()

    assert "row_sample" in df.columns
    group_a_rows = df.filter(pl.col("sample") == "SampleA")
    row_samples = sorted(group_a_rows.get_column("row_sample").unique().to_list())
    assert row_samples == ["SampleA", "SampleA R1", "SampleA R2"]

    reconstructed = ViolinPlotInputData.from_df(df, original.pconfig, original.anchor)
    section_r = list(reconstructed.dt.section_by_id.values())[0]

    assert SampleGroup("SampleA") in section_r.rows_by_sgroup
    group_a = section_r.rows_by_sgroup[SampleGroup("SampleA")]
    assert len(group_a) == 3

    row_sample_names = [str(r.sample) for r in group_a]
    assert "SampleA" in row_sample_names
    assert "SampleA R1" in row_sample_names
    assert "SampleA R2" in row_sample_names

    agg_row = next(r for r in group_a if str(r.sample) == "SampleA")
    assert agg_row.data[ColumnKey("reads")].raw == 20000

    r1_row = next(r for r in group_a if str(r.sample) == "SampleA R1")
    assert r1_row.data[ColumnKey("reads")].raw == 10000

    assert SampleGroup("SampleB") in section_r.rows_by_sgroup
    assert len(section_r.rows_by_sgroup[SampleGroup("SampleB")]) == 1


def test_grouped_samples_survive_merge():
    """
    Verify that merging two ViolinPlotInputData objects with grouped samples
    preserves the group hierarchy and resolved values.
    """
    data1 = _make_grouped_violin_data()

    pconfig = TableConfig(id="test_grouped", title="Test Grouped Table")
    anchor = plot_anchor(pconfig)
    table_anchor = Anchor(f"{anchor}_table")

    grouped_data2: dict = {
        SectionKey("section1"): {
            SampleGroup("SampleC"): [
                InputRow(sample=SampleName("SampleC"), data={"reads": 30000, "quality": 38.0}),
                InputRow(sample=SampleName("SampleC R1"), data={"reads": 15000, "quality": 37.0}),
                InputRow(sample=SampleName("SampleC R2"), data={"reads": 15000, "quality": 39.0}),
            ],
        }
    }
    headers: dict = {
        SectionKey("section1"): {
            ColumnKey("reads"): {"title": "Total Reads"},
            ColumnKey("quality"): {"title": "Quality"},
        }
    }
    dt2 = table_object.DataTable.create(
        data=grouped_data2,
        table_id=pconfig.id,
        table_anchor=table_anchor,
        pconfig=pconfig.model_copy(),
        headers=headers,
    )
    data2 = ViolinPlotInputData(
        dt=dt2,
        plot_type=PlotType.VIOLIN,
        pconfig=pconfig,
        anchor=anchor,
        show_table_by_default=True,
        creation_date=datetime.now() + timedelta(seconds=1),
    )

    merged = ViolinPlotInputData.merge(data1, data2)

    section = list(merged.dt.section_by_id.values())[0]

    assert SampleGroup("SampleA") in section.rows_by_sgroup
    group_a = section.rows_by_sgroup[SampleGroup("SampleA")]
    assert len(group_a) == 3
    agg_a = next(r for r in group_a if str(r.sample) == "SampleA")
    assert agg_a.data[ColumnKey("reads")].raw == 20000

    assert SampleGroup("SampleC") in section.rows_by_sgroup
    group_c = section.rows_by_sgroup[SampleGroup("SampleC")]
    assert len(group_c) == 3
    agg_c = next(r for r in group_c if str(r.sample) == "SampleC")
    assert agg_c.data[ColumnKey("reads")].raw == 30000


def test_grouped_samples_parquet_roundtrip(tmp_path):
    """
    End-to-end test: write grouped sample data to parquet, load it back,
    and verify group hierarchy and resolved values are preserved.
    """
    original = _make_grouped_violin_data()
    original.save_to_parquet()

    flush_to_parquet()
    parquet_path = tmp_dir.parquet_file()
    assert parquet_path.exists()

    df = pl.read_parquet(parquet_path)
    plot_input_rows = df.filter((pl.col("type") == "plot_input") & (pl.col("anchor") == str(original.anchor)))
    assert plot_input_rows.height == 1

    plot_input_json = json.loads(plot_input_rows.get_column("plot_input_data")[0])
    loaded = create_plot_input_data_only(plot_input_json)
    assert isinstance(loaded, ViolinPlotInputData)

    section = list(loaded.dt.section_by_id.values())[0]
    assert SampleGroup("SampleA") in section.rows_by_sgroup
    group_a = section.rows_by_sgroup[SampleGroup("SampleA")]
    assert len(group_a) == 3
    agg = next(r for r in group_a if str(r.sample) == "SampleA")
    assert agg.data[ColumnKey("reads")].raw == 20000
