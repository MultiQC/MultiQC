import pytest

from multiqc.modules.mosdepth.mosdepth import genstats_cov_thresholds, cov_at_cum_fraction, calc_iqr_coverage


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
    if actual_iqr is not None and actual_median is not None: 
        actual_iqr_cv = actual_iqr / actual_median 
        assert actual_iqr_cv == expected_iqr / expected_median
