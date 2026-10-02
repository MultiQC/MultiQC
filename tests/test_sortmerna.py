import multiqc
from multiqc.modules.sortmerna.sortmerna import MultiqcModule


def test_paired_end_uses_first_reads_file(data_dir):
    """Paired-end logs have two 'Reads file' lines, the results belong to the sample named by the first"""
    multiqc.parse_logs(data_dir / "modules/sortmerna/issue_3107")
    mod = MultiqcModule()
    assert list(mod.sortmerna) == ["SRR31139166_chr22_1_val_1"]
    assert mod.sortmerna["SRR31139166_chr22_1_val_1"]["rRNA"] == 1150
