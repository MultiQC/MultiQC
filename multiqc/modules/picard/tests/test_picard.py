import pytest

from multiqc import report
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.picard.picard import MultiqcModule, TOOLS
from multiqc.utils import testing


@pytest.fixture
def data_dir():
    return testing.data_dir()


# Most files are expected to produce one sample, these are exceptions:
NUM_SAMPLES_BY_FILE = {
    "two_samples.txt": 2,  # AlignmentSummaryMetrics
    "single_sample.gvcf.variant_calling_summary_metrics": 0,  # VariantCallingMetrics
    "overview_errors.csv": 0,  # ValidateSamFile
    "overview_warnings.csv": 0,  # ValidateSamFile
    "truncated.txt": 0,  # RnaSeqMetrics
    "P1339_1001.collectInsertSize.txt": 3,  # InsertSizeMetrics
    "zero_reads.InsertSize_metrics.txt": 0,  # InsertSizeMetrics
    "alignment-hs_metrics_2samps.txt": 2,  # HsMetrics
    "A1.sorted.dup.recal.all.metrics.base_distribution_by_cycle_metrics": 2,  # BaseDistributionByCycleMetrics
    "base_distribution_by_cycle.txt": 2,  # BaseDistributionByCycleMetrics
}


@pytest.mark.parametrize("tool", TOOLS)
def test_picard_data_parsed(tool, data_dir):
    tool_subdir = data_dir / "modules/picard" / tool
    print(tool_subdir)
    assert tool_subdir.exists()

    # Check that all test files are parsed
    for path in tool_subdir.rglob("*"):
        if path.name.startswith("."):
            continue
        if path.parent.name == "multiqc_data" or path.name == "multiqc_report.html" or not path.is_file():
            continue

        # print(f"Scanning path {path}")
        report.reset()
        report.analysis_files = [path]
        report.search_files(["picard"])

        try:
            m = MultiqcModule(tools=[tool])
        except ModuleNoSamplesFound:
            samples_parsed = []
        else:
            samples_parsed = list(m.samples_parsed_by_tool.get(tool, []))
        print(f"{path.name}: {samples_parsed}")
        expected_num_samples = NUM_SAMPLES_BY_FILE.get(path.name, 1)
        assert len(samples_parsed) == expected_num_samples
