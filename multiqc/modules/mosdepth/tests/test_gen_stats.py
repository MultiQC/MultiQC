import pytest

from multiqc.modules.mosdepth.mosdepth import genstats_cov_thresholds


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
