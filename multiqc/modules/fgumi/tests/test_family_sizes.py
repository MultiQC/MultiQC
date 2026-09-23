import pytest

from multiqc import report
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.fgumi import MultiqcModule

FAMILY_SIZES_TSV = (
    "family_size\tcount\tfraction\tfraction_gt_or_eq_family_size\n1\t60\t0.6\t1\n2\t30\t0.3\t0.4\n5\t10\t0.1\t0.1\n"
)


def test_family_sizes_parsed_and_plotted(run_fgumi):
    module = run_fgumi({"S1.family_sizes.txt": FAMILY_SIZES_TSV})
    assert module.samples_parsed_by_tool["family_sizes"] == {"S1"}
    assert module.saved_raw_data["multiqc_fgumi_family_sizes"]["S1"] == {1: 60, 2: 30, 5: 10}
    assert "fgumi_family_sizes" in report.plot_by_id


def test_user_named_histogram_is_detected_by_header(run_fgumi):
    module = run_fgumi({"my_hist.tsv": FAMILY_SIZES_TSV})
    assert module.samples_parsed_by_tool["family_sizes"] == {"my_hist"}


def test_sample_name_strips_fgumi_suffixes(run_fgumi):
    # A second file with a longer fgumi suffix must map to the same sample, not "S1.duplex_family_sizes".
    from multiqc.modules.fgumi.util import sample_name

    module = run_fgumi({"S1.family_sizes.txt": FAMILY_SIZES_TSV})
    for fn in ["S1.duplex_family_sizes.txt", "S1.grouping_metrics.txt", "S1.umi_counts.txt.gz", "S1.family_sizes.txt"]:
        f = {"fn": fn, "root": "", "s_name": fn, "sp_key": "fgumi/family_sizes"}  # the keys find_log_files provides
        assert sample_name(module, f) == "S1", fn


def test_family_sizes_claimed_by_fgumi_not_fgbio(tmp_path):
    from multiqc.modules.fgbio import MultiqcModule as FgbioModule

    path = tmp_path / "S1.family_sizes.txt"
    path.write_text(FAMILY_SIZES_TSV)
    report.reset()
    report.analysis_files = [path]
    report.search_files(["fgbio", "fgumi"])
    with pytest.raises(ModuleNoSamplesFound):
        FgbioModule()
    assert MultiqcModule().samples_parsed_by_tool["family_sizes"] == {"S1"}


def test_malformed_file_is_skipped(run_fgumi, caplog):
    truncated = FAMILY_SIZES_TSV + "7\t3\n"
    with pytest.raises(ModuleNoSamplesFound):
        run_fgumi({"S1.family_sizes.txt": truncated})
    assert any("S1.family_sizes.txt" in record.message for record in caplog.records)


def test_no_fgumi_files_raises(run_fgumi):
    with pytest.raises(ModuleNoSamplesFound):
        run_fgumi({"irrelevant.txt": "hello world\n"})
