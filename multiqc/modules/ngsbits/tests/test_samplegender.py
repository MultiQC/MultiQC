from multiqc import report
from multiqc.modules.ngsbits import MultiqcModule


def test_parse(tmp_path):
    f1 = tmp_path / "sample1_ngsbits_sex.tsv"
    f1.write_text("""\
#file   gender  reads_chry      reads_chrx      ratio_chry_chrx
sample_sorted_md.bam        male    3872086 14244970        0.2718
""")

    report.analysis_files = [f1]
    report.search_files(["*ngsbits_sex.tsv"])
    m = MultiqcModule()
    assert len(m.saved_raw_data) == 1
    print(m.saved_raw_data)

    assert len(m.sections) == 1
    assert len(report.general_stats_data) > 0
