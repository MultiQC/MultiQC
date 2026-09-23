from multiqc import report
from multiqc.modules.fgumi.umis import summarize_observations

UMI_COUNTS = (
    "umi\traw_observations\traw_observations_with_errors\tunique_observations\tfraction_raw_observations"
    "\tfraction_unique_observations\n"
    "AAAA\t1\t0\t1\t0.1\t0.25\nCCCC\t1\t0\t1\t0.1\t0.25\nGGGG\t3\t1\t1\t0.3\t0.25\nTTTT\t5\t0\t1\t0.5\tNaN\n"
)
CORRECT = (
    "umi\ttotal_matches\tperfect_matches\tone_mismatch_matches\ttwo_mismatch_matches\tother_matches"
    "\tfraction_of_matches\trepresentation\n"
    "ACGT\t90\t80\t8\t2\t0\t0.9\t1.8\nTTTT\t0\t0\t0\t0\t0\t0\t0\nNNNN\t10\t0\t0\t0\t10\t0.1\tInfinity\n"
)


def test_summarize_observations_median_singletons_and_log2_bins():
    summary = summarize_observations([1, 1, 3, 5])
    assert summary["n"] == 4
    assert summary["median"] == 2.0
    assert summary["singleton_pct"] == 50.0
    assert summary["bins"] == {1: 2, 4: 1, 8: 1}  # bins are (1], (2], (3-4], (5-8], ... keyed by upper bound


def test_umi_counts_summary(run_fgumi):
    module = run_fgumi({"S1.umi_counts.txt": UMI_COUNTS})
    assert module.saved_raw_data["multiqc_fgumi_umi_counts"]["S1"] == {"n": 4, "median": 2.0, "singleton_pct": 50.0}
    assert "fgumi_umi_counts" in report.plot_by_id


def test_umi_correction_summary_counts_unmatched_all_n_row(run_fgumi):
    module = run_fgumi({"S1.txt": CORRECT})
    assert module.saved_raw_data["multiqc_fgumi_correct"]["S1"] == {
        "perfect": 80,
        "one_mismatch": 8,
        "two_mismatch": 2,
        "other": 0,
        "unmatched": 10,
    }


DUPLEX_UMI_COUNTS = (
    "umi\traw_observations\traw_observations_with_errors\tunique_observations\tfraction_raw_observations"
    "\tfraction_unique_observations\tfraction_unique_observations_expected\n"
    "AAAA-CCCC\t2\t0\t1\t0.5\t0.5\t0.4\nCCCC-AAAA\t2\t0\t1\t0.5\t0.5\t0.4\n"
)


def test_duplex_umi_counts_summary(run_fgumi):
    module = run_fgumi({"S1.duplex_umi_counts.txt": DUPLEX_UMI_COUNTS})
    assert module.saved_raw_data["multiqc_fgumi_duplex_umi_counts"]["S1"] == {
        "n": 2,
        "median": 2.0,
        "singleton_pct": 0.0,
    }
    assert "fgumi_duplex_umi_counts" in report.plot_by_id
