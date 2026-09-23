import pytest

from multiqc import report

HEADER = (
    "sample\tlibrary\tfiltered_templates\tfiltered_malformed_record\tfiltered_no_primary_reads\tfiltered_unmapped"
    "\tfiltered_not_passing_filter\tfiltered_low_mapping_quality\tfiltered_low_mate_mapping_quality"
    "\tfiltered_missing_umi\tfiltered_ns_in_umi\tfiltered_umi_too_short\tpassthrough_templates\ttotal_templates"
    "\tunique_templates\tduplicate_templates\tduplicate_rate\ttotal_reads\tunique_reads\tduplicate_reads"
    "\tsecondary_reads\tsupplementary_reads\tmissing_tc_tag\tmapped_pairs\tduplicate_pairs\tmapped_orphans"
    "\tduplicate_orphans\tunmapped_pairs\tunmapped_orphans\tunmated_templates\tpercent_duplication"
    "\testimated_library_size\n"
)
ZEROS = "\t".join(["0"] * 11)  # filtered_templates .. passthrough_templates
DEDUP = (
    HEADER
    + f"S1\tlibA\t{ZEROS}\t100\t80\t20\t0.2\t200\t160\t40\t0\t0\t0\t100\t20\t0\t0\t0\t0\t0\t0.2\t500\n"
    + f"S1\tAll Reads\t{ZEROS}\t100\t80\t20\t0.2\t200\t160\t40\t0\t0\t0\t100\t20\t0\t0\t0\t0\t0\tNaN\t\n"
)
LADDER = (
    "library\ttemplates_seen\tduplicate_fraction\twindow_templates\twindow_duplicate_fraction\n"
    "libA\t50\t0.1\t50\t0.1\nlibA\t100\t0.2\t50\t0.3\n"
)


def test_dedup_general_stats_from_all_reads_row_and_non_finite(run_fgumi):
    module = run_fgumi({"S1.txt": DEDUP})
    assert module.samples_parsed_by_tool["dedup"] == {"S1"}
    assert module.saved_raw_data["multiqc_fgumi_dedup"]["S1"] == {
        "unique_templates": 80,
        "duplicate_templates": 20,
        "percent_duplication": None,  # NaN must become None, not 0
        "estimated_library_size": None,  # empty cell
    }
    assert "fgumi_dedup_templates" in report.plot_by_id


def test_duplication_ladder(run_fgumi):
    module = run_fgumi({"S1.txt": LADDER})
    ladder = module.saved_raw_data["multiqc_fgumi_dedup_ladder"]["S1"]
    assert ladder["cumulative"] == pytest.approx({50: 10.0, 100: 20.0})
    assert ladder["window"] == pytest.approx({50: 10.0, 100: 30.0})


TWO_LIBRARY_LADDER = LADDER + "libB\t50\t0.2\t50\t0.2\n"


def test_ladder_honours_sample_filters_on_the_sample_name(run_fgumi):
    from multiqc import config

    original = config.sample_names_ignore
    config.sample_names_ignore = ["S1"]
    try:
        module = run_fgumi({"S1.txt": TWO_LIBRARY_LADDER, "S2.txt": LADDER})
    finally:
        config.sample_names_ignore = original
    assert set(module.saved_raw_data["multiqc_fgumi_dedup_ladder"]) == {"S2"}


def test_dedup_without_all_reads_row_warns(run_fgumi, caplog):
    only_library = "".join(line + "\n" for line in DEDUP.splitlines() if "All Reads" not in line)
    run_fgumi({"S1.txt": only_library, "S2.txt": LADDER})
    assert any("S1.txt" in r.message and "All Reads" in r.message for r in caplog.records)
