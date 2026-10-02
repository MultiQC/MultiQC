from multiqc import report
from multiqc.modules.checkm2 import MultiqcModule

ALLMODELS_REPORT = (
    "Name\tCompleteness_General\tContamination\tCompleteness_Specific\tCompleteness_Model_Used\tAdditional_Notes\n"
    "test1.faa\t26.36\t0.0\t15.9\tNeural Network (Specific Model)\tNone\n"
    "test2.faa\t96.06\t0.05\t99.44\tNeural Network (Specific Model)\tNone\n"
)


def test_allmodels_report_is_detected(tmp_path):
    (tmp_path / "checkm2-report.tsv").write_text(ALLMODELS_REPORT)

    report.reset()
    report.analysis_files = [tmp_path]
    report.search_files(["checkm2"])
    MultiqcModule()

    samples = {str(s) for section in report.general_stats_data.values() for s in section}
    assert samples == {"test1", "test2"}
    assert any("Completeness_General" in str(h) for h in next(iter(report.general_stats_headers.values())))
