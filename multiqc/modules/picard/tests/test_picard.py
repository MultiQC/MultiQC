import pytest

from multiqc import report
from multiqc.modules.picard.picard import MultiqcModule, TOOLS
from multiqc.utils import testing


@pytest.fixture
def data_dir():
    return testing.data_dir()


# Most files are expected to produce one sample, these are exceptions:
NUM_SAMPLES_BY_FILE = {
    "two_samples.txt": 2,
}


@pytest.mark.parametrize("tool", TOOLS)
def test_picard_data_parsed(tool, data_dir):
    tool_subdir = data_dir / "modules/picard" / tool
    print(tool_subdir)
    assert tool_subdir.exists()

    # Check that all test files are parsed
    for path in tool_subdir.rglob("*"):
        if path.parent.name == "multiqc_data" or path.name == "multiqc_report.html" or not path.is_file():
            continue

        print(f"Scanning path {path}")
        report.reset()
        report.analysis_files = [path]
        report.search_files(["picard"])
        m = MultiqcModule(tools=[tool])
        expected_num_samples = NUM_SAMPLES_BY_FILE.get(path.name, 1)
        assert len(m.saved_raw_data) == 1
        samples_parsed = list(m.saved_raw_data.values())[0].keys()
        print(path, samples_parsed)
        assert len(samples_parsed) == expected_num_samples
