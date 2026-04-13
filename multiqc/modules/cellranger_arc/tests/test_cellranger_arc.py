import glob

import pytest
from multiqc import config, report
from multiqc.modules.cellranger_arc import MultiqcModule
from multiqc.utils import testing


@pytest.fixture
def data_dir():
    return testing.data_dir()


def test_cellranger_arc(data_dir):
    html_files = glob.glob(str(data_dir / "modules/cellranger_arc/arc-*/*.html"))
    for file in html_files:
        report.analysis_files = [file]
        report.search_files(["cellranger_arc"])

    config.preserve_module_raw_data = True
    m = MultiqcModule()
    assert m.saved_raw_data is not None
    assert len(m.saved_raw_data) > 0

    ## check if both sameples are present in the data
    data = m.saved_raw_data["multiqc_cellranger_arc"]
    assert "TEST_SAMPLE" in data, "Missing sample: 'TEST_SAMPLE'"
    assert "P210059__KIZH_03" in data, "Missing sample: 'P210059__KIZH_03'"

    ## check if general stats data is there and there are 8 sections
    assert len(report.general_stats_data) > 0
    assert len(m.sections) == 8

    ## check of the plots have data
    for sam in ["TEST_SAMPLE", "P210059__KIZH_03"]:
        for key in m.plots_data_by_sample:
            assert m.plots_data_by_sample[key][sam], f"{key} plot data is empty for {sam}"
