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


@pytest.mark.parametrize(
    "fn",
    [
        "S1.family_sizes.txt",
        "S1.duplex_family_sizes.txt",  # a longer suffix must not leave "S1.duplex_family_sizes"
        "S1.grouping_metrics.txt",
        "S1.umi_counts.txt",
        # `fgumi runall --all-metrics S1` puts the stage before the suffix
        "S1.group.family_sizes.txt",
        "S1.duplex.umi_counts.txt",
        "S1.simplex.simplex_yield_metrics.txt",
        "S1.codec.family_sizes.txt",
        "S1.correct.metrics.txt",
        "S1.filter.stats.txt",
    ],
)
def test_sample_name_strips_fgumi_suffixes(run_fgumi, fn):
    from multiqc.modules.fgumi.util import sample_name

    module = run_fgumi({"S1.family_sizes.txt": FAMILY_SIZES_TSV})
    f = {"fn": fn, "root": "", "s_name": fn, "sp_key": "fgumi/family_sizes"}  # the keys find_log_files provides
    assert sample_name(module, f) == "S1"


def test_fullnames_keeps_the_file_name(run_fgumi):
    from multiqc import config

    config.fn_clean_sample_names = False
    module = run_fgumi({"S1.family_sizes.txt": FAMILY_SIZES_TSV})
    assert module.samples_parsed_by_tool["family_sizes"] == {"S1.family_sizes.txt"}


def test_runall_outputs_share_one_sample(run_fgumi):
    grouping = "\t".join(
        ["accepted_sam_records", "discarded_non_pf", "discarded_poor_alignment", "discarded_ns_in_umi"]
        + ["discarded_umis_to_short"]
    )
    module = run_fgumi(
        {
            "S1.group.family_sizes.txt": FAMILY_SIZES_TSV,
            "S1.group.grouping_metrics.txt": grouping + "\n10\t0\t0\t0\t0\n",
        }
    )
    assert module.samples_parsed_by_tool["family_sizes"] == {"S1"}
    assert module.samples_parsed_by_tool["grouping"] == {"S1"}
    # Each file is listed as a source for the sample; neither replaces the other.
    sources = {s for s, by_sample in report.data_sources["fgumi"].items() if "S1" in by_sample}
    assert {"family_sizes", "grouping_metrics"} <= sources


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
