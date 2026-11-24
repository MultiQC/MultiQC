import pytest

from multiqc import report
from multiqc.modules.mosdepth import MultiqcModule
from multiqc.modules.mosdepth.mosdepth import calc_iqr_coverage, cov_at_cum_fraction, genstats_cov_thresholds
from multiqc.types import SampleGroup


def test_genstats_cov_thresholds():
    cum_fraction_by_cov = {
        1: 1.0,
        10: 0.8,
        20: 0.2,
        30: 0.1,
    }
    thresholds = 10, 15, 30, 200

    actual_thresholds = genstats_cov_thresholds(cum_fraction_by_cov, thresholds)
    assert actual_thresholds == {
        "10_x_pc": 80.0,
        "15_x_pc": 20.0,
        "30_x_pc": 10.0,
        "200_x_pc": 0.0,
    }


@pytest.mark.parametrize(
    "cum_fraction_by_cov,expected_median,expected_iqr",
    [
        ({1: 1.0, 10: 0.8, 20: 0.3, 30: 0.1}, 10, 10),
        ({}, None, None),
        ({1: 1.0}, 1, 0),
    ],
)
def test_calc_median_and_iqr_coverage(cum_fraction_by_cov, expected_median, expected_iqr):
    actual_median = cov_at_cum_fraction(cum_fraction_by_cov, 0.5)
    assert actual_median == expected_median
    actual_iqr = calc_iqr_coverage(cum_fraction_by_cov)
    assert actual_iqr == expected_iqr


def test_iqr_cv_in_general_stats(tmp_path):
    dist_file = tmp_path / "sampleA.mosdepth.global.dist.txt"
    dist_file.write_text("total\t30\t0.10\ntotal\t20\t0.30\ntotal\t10\t0.80\ntotal\t1\t1.00\n")

    report.reset()
    report.analysis_files = [dist_file]
    report.search_files(["mosdepth"])
    MultiqcModule()

    headers = next(iter(report.general_stats_headers.values()))
    assert "coverage_iqr_cv" in headers
    rows = next(iter(report.general_stats_data.values()))
    # median 10, Q1 10, Q3 20
    assert rows[SampleGroup("sampleA")][0].data["coverage_iqr_cv"] == 1.0
