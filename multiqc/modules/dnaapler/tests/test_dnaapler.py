from multiqc import report
from multiqc.modules.dnaapler import MultiqcModule


def test_dnaapler(data_dir):
    report.reset()
    report.analysis_files = [data_dir / "modules/dnaapler"]
    report.search_files(["dnaapler"])

    module = MultiqcModule()
    assert len(module.sections) == 2
    assert len(module.general_stats_data) > 0
