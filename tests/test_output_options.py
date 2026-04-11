import json
import os

import polars as pl
import pytest

import multiqc
from multiqc import config, report, BaseMultiqcModule, write_report
from multiqc.core.update_config import ClConfig
from multiqc.plots import table
from multiqc.types import Anchor


@pytest.fixture()
def stub_modules():
    """
    Set stub modules to make write_report work
    """
    report.modules = [BaseMultiqcModule()]


@pytest.mark.parametrize(
    "options,expected_files",
    [
        ({}, {"multiqc_report.html", "multiqc_data"}),
        ({"filename": "NAME"}, {"NAME.html", "NAME_data"}),
        ({"filename": "NAME.html"}, {"NAME.html", "NAME_data"}),
        ({"title": "My Title"}, {"My-Title_multiqc_report.html", "My-Title_multiqc_report_data"}),
        ({"make_data_dir": False}, {"multiqc_report.html"}),
        ({"make_report": False}, {"multiqc_data"}),
        ({"make_report": False, "make_data_dir": False}, set()),
    ],
)
def test_filename(stub_modules, tmp_path, options, expected_files):
    """
    Verify that the filename option works
    """
    write_report(output_dir=tmp_path, **options)

    assert set(os.listdir(tmp_path)) == expected_files


def test_write_stdout(stub_modules, tmp_path, capsys):
    """
    Verify stdout option: stdout contains only HTML, nothing else is written to disk
    """
    files_before = set(os.listdir(tmp_path))
    os.chdir(tmp_path)

    write_report(filename="stdout")

    captured = capsys.readouterr()
    assert captured.out.startswith("""<!doctype html>\n<html lang="en">\n  <head>""")

    files_after = set(os.listdir(tmp_path))
    assert files_before == files_after


def test_zip_data_dir_with_rerun(stub_modules, tmp_path):
    """
    Verify that zip_data_dir works correctly when re-running with --force.
    This tests the fix for the bug where a file named 'multiqc_data' was created,
    causing NotADirectoryError on subsequent runs.
    See https://github.com/MultiQC/MultiQC/issues/3358
    """
    # First run with zip_data_dir
    write_report(output_dir=tmp_path, zip_data_dir=True)

    # Check that the zip file was created and the data directory was removed
    assert (tmp_path / "multiqc_data.zip").exists()
    assert not (tmp_path / "multiqc_data").exists()

    # Second run with --force should work without errors
    write_report(output_dir=tmp_path, zip_data_dir=True, force=True)

    # Check again
    assert (tmp_path / "multiqc_data.zip").exists()
    assert not (tmp_path / "multiqc_data").exists()


def test_parquet_wide_merges_samples(tmp_path):
    """
    Verify that wide parquet format merges multiple tables for the same sample into a single row.
    Verify that the values all get persisted in the proper format.
    """
    with open(tmp_path / "f1_mqc.json", "w") as f:
        json.dump(
            {
                "data": {"foo": {"col1": None, "col2": 3.4, "col3": float("nan"), "col4": "bar", "col5": True}},
                "id": "myid",
                "anchor": "myanchor",
                "plot_type": "table",
                "file_format": "json",
                "section_name": "My Section",
                "description": "my desc",
            },
            f,
        )
    with open(tmp_path / "f2_mqc.json", "w") as f:
        json.dump(
            {
                "data": {"foo": {"col1": 2}},
                "id": "myid2",
                "anchor": "myanchor2",
                "plot_type": "table",
                "file_format": "json",
                "section_name": "My Section",
                "description": "my desc",
            },
            f,
        )
    with open(tmp_path / "f3_mqc.json", "w") as f:
        json.dump(
            {
                "data": {"foo": {"col1": 3}},
                "id": "myid3",
                "anchor": "myanchor3",
                "plot_type": "table",
                "file_format": "json",
                "section_name": "My Section",
                "description": "my desc",
            },
            f,
        )

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    multiqc.run(tmp_path, cfg=ClConfig(output_dir=out_dir, cl_config=["parquet_format: wide"]))

    df = pl.read_parquet(out_dir / "multiqc_data" / "multiqc.parquet")

    table_rows = df.filter(pl.col("type") == "table_row")

    # All three tables should be merged into a single row for sample "foo"
    assert table_rows.height == 1
    assert (table_rows["sample"] == "foo").all()

    # Each table should have contributed its metric column
    assert (table_rows["myid / col1"].is_null()).all()
    assert (table_rows["myid2 / col1"] == 2).all()
    assert (table_rows["myid3 / col1"] == 3).all()

    # typing tests
    assert (table_rows["myid / col2"] == 3.4).all()
    assert (table_rows["myid / col3"].is_nan()).all()
    assert (table_rows["myid / col4"] == "bar").all()
    assert (table_rows["myid / col5"] == True).all()


def test_parquet_wide_persists_modified_values(tmp_path):
    """
    Verify that the post-modify value is what gets persisted in wide parquet, not the raw value.
    """
    config.parquet_format = "wide"
    module = BaseMultiqcModule(name="test-module", anchor=Anchor("test_module"))
    module.add_section(
        name="Test Section",
        plot=table.plot(
            data={"sample1": {"metric1": 5}},
            headers={"metric1": {"title": "Metric 1", "modify": lambda x: x * 10}},
            pconfig={"id": "test_table", "title": "Test: Test Table"},
        ),
    )
    report.modules = [module]
    write_report(force=True, output_dir=str(tmp_path))

    df = pl.read_parquet(tmp_path / "multiqc_data" / "multiqc.parquet")
    table_rows = df.filter(pl.col("type") == "table_row")
    assert table_rows.height == 1
    assert (table_rows["sample"] == "sample1").all()
    assert (table_rows["test_table / metric1"] == 50.0).all()
