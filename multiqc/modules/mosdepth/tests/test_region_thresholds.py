import gzip
import logging

import pytest

from multiqc import config, report
from multiqc.modules.mosdepth import MultiqcModule
from multiqc.modules.mosdepth.mosdepth import count_passing_regions

HEADER = "#chrom\tstart\tend\tregion\t1X\t10X\t20X\n"


@pytest.fixture(autouse=True)
def clean_report():
    report.reset()
    yield
    report.reset()


def test_count_passing_regions():
    lines = [
        HEADER,
        "1\t0\t100\tA\t100\t100\t100\n",  # passes everywhere
        "1\t100\t200\tA\t100\t99\t50\n",  # exactly 99% at 10X passes, 50% at 20X fails
        "1\t200\t300\tB\t0\t0\t0\n",  # fails everywhere
        "1\t300\t300\tC\t0\t0\t0\n",  # zero-length, left out
    ]
    n_regions, passing = count_passing_regions(lines, 99.0)
    assert n_regions == 3
    assert passing == {1: 2, 10: 2, 20: 1}


def test_count_passing_regions_custom_pct():
    lines = [HEADER, "1\t0\t100\tA\t100\t90\t50\n"]
    assert count_passing_regions(lines, 90.0) == (1, {1: 1, 10: 1, 20: 0})


@pytest.mark.parametrize(
    "lines",
    [
        [],
        ["#chrom\tstart\tend\t1X\n"],
        [HEADER, "1\t0\t100\tA\t100\t100\n"],
        [HEADER, "1\t0\t100\tA\t100\tx\t100\n"],
    ],
    ids=["empty", "bad_header", "short_row", "non_numeric"],
)
def test_count_passing_regions_rejects_bad_input(lines):
    with pytest.raises(ValueError):
        count_passing_regions(lines, 99.0)


def write_gz(path, text):
    with gzip.open(path, "wt") as fh:
        fh.write(text)


def run_module(tmp_path):
    report.analysis_files = [tmp_path]
    report.search_files(["mosdepth"])
    MultiqcModule()
    return {
        str(sample): values
        for section in report.general_stats_data.values()
        for sample, rows in section.items()
        for values in [rows[0].data]
    }


def test_general_stats_columns(tmp_path):
    write_gz(
        tmp_path / "good.thresholds.bed.gz",
        HEADER + "1\t0\t100\tA\t100\t100\t100\n1\t100\t200\tB\t100\t100\t0\n",
    )
    stats = run_module(tmp_path)
    assert set(stats) == {"good"}
    good = stats["good"]
    assert good["regions_total"] == 2
    assert good["10_x_regions_pc"] == 100.0
    assert good["10_x_regions_failed"] == 0
    assert good["20_x_regions_pc"] == 50.0
    assert good["20_x_regions_failed"] == 1


def test_region_pass_pct_config(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "mosdepth_config", {"region_pass_pct": 50}, raising=False)
    write_gz(tmp_path / "good.thresholds.bed.gz", HEADER + "1\t0\t100\tA\t100\t60\t40\n")
    good = run_module(tmp_path)["good"]
    assert good["10_x_regions_failed"] == 0
    assert good["20_x_regions_failed"] == 1


def test_zero_threshold_column_is_dropped(tmp_path):
    write_gz(tmp_path / "good.thresholds.bed.gz", "#chrom\tstart\tend\tregion\t0X\t1X\n1\t0\t100\tA\t100\t100\n")
    good = run_module(tmp_path)["good"]
    assert "0_x_regions_pc" not in good
    assert good["1_x_regions_pc"] == 100.0


@pytest.mark.parametrize(
    "write_bad",
    [
        lambda p: p.write_text(HEADER),  # not gzipped
        lambda p: write_gz(p, "chr1\t0\t100\tpeak1\t5.0\n"),  # another tool's BED
        lambda p: p.write_bytes(gzip.compress(HEADER.encode())[:-10]),  # truncated archive
    ],
    ids=["not_gzip", "other_tool", "truncated"],
)
def test_unusable_file_is_skipped(tmp_path, caplog, write_bad):
    write_gz(tmp_path / "good.thresholds.bed.gz", HEADER + "1\t0\t100\tA\t100\t100\t100\n")
    bad = tmp_path / "bad.thresholds.bed.gz"
    write_bad(bad)
    with caplog.at_level(logging.WARNING):
        stats = run_module(tmp_path)
    assert f"Skipping {bad}" in caplog.text
    assert set(stats) == {"good"}
