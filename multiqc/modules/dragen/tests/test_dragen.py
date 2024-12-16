import pytest

from multiqc import report
from multiqc.modules.dragen.dragen import MultiqcModule
from multiqc.utils import testing


@pytest.fixture
def data_dir():
    return testing.data_dir()


EXPECTED_SAMPLES_BY_TOOL = {
    "add_mapping_metrics": 5,
    "add_vc_metrics": 3,
    "add_ploidy_estimation_metrics": 2,
    "collect_overall_mean_cov_data": 5,
    "add_coverage_metrics": 15,
    "add_coverage_hist": 3,
    "add_coverage_per_contig": 3,
    "add_fragment_length_hist": 3,
    "add_gc_metrics_hist": 1,
    "add_trimmer_metrics": 1,
    "add_time_metrics": 2,
    "add_rna_metrics": 0,
    "add_rna_transcript_coverage": 0,
    "add_sc_rna_metrics": 0,
    "add_sc_atac_metrics": 0,
}


def test_dragen_data_parsed(data_dir):
    report.reset()
    report.analysis_files = [data_dir / "modules/dragen"]
    report.search_files(["dragen"])

    m = MultiqcModule()
    for tool, samples in m.samples_parsed_by_tool.items():
        assert tool in EXPECTED_SAMPLES_BY_TOOL
        assert len(samples) == EXPECTED_SAMPLES_BY_TOOL[tool]
