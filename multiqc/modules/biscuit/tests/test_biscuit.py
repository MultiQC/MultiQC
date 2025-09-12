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

    parsed = parse_align_mapq(mapq_table)

    # Four keys are expected (frac_align, opt_align, sub_align, mapqs)
    assert len(parsed) == 4

    # mapqs key should have 61 keys
    assert len(parsed["mapqs"]) == 61

    # Check correctly parsed values
    assert parsed["frac_align"] - 98.81 < EPSILON
    assert parsed["opt_align"] == 56535091
    assert parsed["sub_align"] == 9287487
    assert parsed["not_align"] == 794210

    # Spot check MAPQ values
    assert parsed["mapqs"][0] - 9.41 < EPSILON
    assert parsed["mapqs"][14] - 0.05 < EPSILON
    assert parsed["mapqs"][26] - 0.04 < EPSILON
    assert parsed["mapqs"][40] - 3.83 < EPSILON
    assert parsed["mapqs"][60] - 78.25 < EPSILON
