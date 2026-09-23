import pytest

from multiqc import report
from multiqc.base_module import ModuleNoSamplesFound

STATS = (
    "key\tvalue\tdescription\n"
    "raw_reads_considered\t1000\tTotal raw reads considered\n"
    "raw_reads_rejected\t200\tRejected reads\n"
    "raw_reads_used\t800\tUsed reads\n"
    "frac_raw_reads_used\t0.8\tFraction used\n"
    "raw_reads_rejected_for_insufficient_support\t150\tToo few reads\n"
    "raw_reads_rejected_for_minority_alignment\t50\tMinority alignment\n"
    "raw_reads_rejected_for_orphan_consensus\t0\tOrphans\n"
    "consensus_reads_emitted\t300\tConsensus reads\n"
)


def test_consensus_stats(run_fgumi):
    module = run_fgumi({"S1.txt": STATS})
    assert module.samples_parsed_by_tool["consensus_stats"] == {"S1"}
    saved = module.saved_raw_data["multiqc_fgumi_consensus_stats"]["S1"]
    assert {
        k: saved[k]
        for k in (
            "raw_reads_used",
            "raw_reads_rejected_for_insufficient_support",
            "raw_reads_rejected_for_minority_alignment",
        )
    } == {
        "raw_reads_used": 800.0,
        "raw_reads_rejected_for_insufficient_support": 150.0,
        "raw_reads_rejected_for_minority_alignment": 50.0,
    }
    assert "fgumi_consensus_rejections" in report.plot_by_id


def test_non_fgumi_kv_file_is_ignored(run_fgumi):
    with pytest.raises(ModuleNoSamplesFound):
        run_fgumi({"other.tsv": "key\tvalue\tdescription\ncolor\tblue\ta color\n"})


def test_consensus_stats_data_file_keeps_every_parsed_value(run_fgumi):
    module = run_fgumi({"S1.txt": STATS})
    saved = module.saved_raw_data["multiqc_fgumi_consensus_stats"]["S1"]
    assert saved["raw_reads_considered"] == 1000.0
    assert saved["consensus_reads_emitted"] == 300.0
    assert saved["raw_reads_rejected_for_orphan_consensus"] == 0.0
