import pytest

from multiqc import config, report
from multiqc.modules.biscuit import MultiqcModule
from multiqc.utils import testing

# NOTE: These tests could be fleshed out more (more inputs, more error cases, etc.)

EPSILON = 0.0001


@pytest.fixture
def data_dir():
    return testing.data_dir()


def read_file(data_dir, fname):
    with (data_dir / "modules/biscuit/v0.3.16.20200420" / fname).open() as fh:
        return fh.read()


@pytest.fixture
def mapq_table(data_dir):
    return read_file(data_dir, "tcga_lusc_normal_subsampled_mapq_table.txt")


def test_mapq_table(mapq_table):
    """Test that parsing the MAPQ table works as expected"""
    from multiqc.modules.biscuit.biscuit import parse_align_mapq

    parsed = parse_align_mapq(mapq_table, "tcga_lusc_normal_subsampled_mapq_table.txt")

    # Four keys are expected (frac_align, opt_align, sub_align, not_align, mapqs)
    assert len(parsed) == 5

    # mapqs key should have 61 keys
    assert len(parsed["mapqs"]) == 61

    # Check correctly parsed values
    assert abs(parsed["frac_align"] - 98.80779) < EPSILON
    assert parsed["opt_align"] == 56535091
    assert parsed["sub_align"] == 9287487
    assert parsed["not_align"] == 794210

    # Spot check MAPQ values
    assert abs(parsed["mapqs"][0] - 9.41347) < EPSILON
    assert abs(parsed["mapqs"][14] - 0.05310) < EPSILON
    assert abs(parsed["mapqs"][26] - 0.04399) < EPSILON
    assert abs(parsed["mapqs"][40] - 3.82724) < EPSILON
    assert abs(parsed["mapqs"][60] - 78.25065) < EPSILON
