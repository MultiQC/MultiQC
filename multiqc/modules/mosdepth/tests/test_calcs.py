import pytest

from multiqc.modules.mosdepth.mosdepth import genstats_cov_thresholds, calc_median_coverage, calc_iqr_coverage


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


def test_calc_median_and_iqr_coverage():
    cum_fraction_by_cov = {
        1: 1.0,
        10: 0.8,
        20: 0.3,
        30: 0.1,
    }

    actual_median = calc_median_coverage(cum_fraction_by_cov)
    assert actual_median == 10
    actual_iqr = calc_iqr_coverage(cum_fraction_by_cov)
    assert actual_iqr == 10
    actual_iqr_cv = actual_iqr / actual_median
    assert actual_iqr_cv == 1.0
