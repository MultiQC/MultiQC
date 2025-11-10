import pytest

from multiqc import config, report
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.ribotish import MultiqcModule
from multiqc.utils import testing


@pytest.fixture
def data_dir():
    return testing.data_dir()


def test_ribotish_data_parsed(data_dir):
    """Test that RiboTish qual files are parsed correctly"""
    data_subdir = data_dir / "modules/ribotish"

    # Find all qual.txt files
    qual_files = list(data_subdir.glob("*_qual.txt"))
    assert len(qual_files) > 0, "No qual.txt files found in test data"

    # Test parsing each file
    for path in qual_files:
        report.reset()
        report.analysis_files = [path]
        report.search_files(["ribotish"])
        config.preserve_module_raw_data = True

        try:
            m = MultiqcModule()
        except ModuleNoSamplesFound:
            pytest.fail(f"Failed to parse {path.name}")

        # Check that data was parsed
        assert m.ribotish_data is not None
        assert len(m.ribotish_data) == 1, f"Expected 1 sample from {path.name}, got {len(m.ribotish_data)}"

        # Check that frame proportions were calculated
        assert m.frame_proportions is not None
        assert len(m.frame_proportions) == 1


def test_ribotish_frame_calculation(data_dir):
    """Test that frame proportion calculations are correct"""
    data_subdir = data_dir / "modules/ribotish"
    qual_file = next(data_subdir.glob("*_qual.txt"))

    report.reset()
    report.analysis_files = [qual_file]
    report.search_files(["ribotish"])

    m = MultiqcModule()

    # Get the sample name (should be only one)
    sample_name = list(m.frame_proportions.keys())[0]

    # Check that frame proportions sum to 1.0 for each read length
    for length, props in m.frame_proportions[sample_name].items():
        total_prop = props["f0_prop"] + props["f1_prop"] + props["f2_prop"]
        assert abs(total_prop - 1.0) < 0.001 or props["total"] == 0, (
            f"Frame proportions don't sum to 1.0 for length {length}: {total_prop}"
        )


def test_ribotish_multiple_samples(data_dir):
    """Test parsing multiple samples at once"""
    data_subdir = data_dir / "modules/ribotish"
    qual_files = list(data_subdir.glob("*_qual.txt"))

    report.reset()
    report.analysis_files = qual_files
    report.search_files(["ribotish"])

    m = MultiqcModule()

    # Should have parsed all samples
    assert len(m.ribotish_data) == len(qual_files), f"Expected {len(qual_files)} samples, got {len(m.ribotish_data)}"

    # Check that all samples have frame proportions calculated
    assert len(m.frame_proportions) == len(qual_files)
