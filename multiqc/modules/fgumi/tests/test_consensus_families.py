from multiqc import report

SIMPLEX_FS = (
    "family_size\tcs_count\tcs_fraction\tcs_fraction_gt_or_eq_size\tss_count\tss_fraction\tss_fraction_gt_or_eq_size\n"
    "1\t50\t0.5\t1\t40\t0.4\t1\n"
    "3\t50\t0.5\t0.5\t60\t0.6\tNaN\n"
)
DUPLEX_FS = (
    "family_size\tcs_count\tcs_fraction\tcs_fraction_gt_or_eq_size\tss_count\tss_fraction\tss_fraction_gt_or_eq_size"
    "\tds_count\tds_fraction\tds_fraction_gt_or_eq_size\n"
    "1\t50\t0.5\t1\t40\t0.4\t1\t30\t0.3\t1\n"
    "2\t50\t0.5\t0.5\t60\t0.6\t0.6\t70\t0.7\tInfinity\n"
)
DUPLEX_AB_BA = (
    "ab_size\tba_size\tcount\tfraction\tfraction_gt_or_eq_size\n"
    "1\t0\t40\t0.4\t1\n"
    "2\t1\t30\t0.3\t0.6\n"
    "25\t22\t30\t0.3\t0.3\n"
)


def test_simplex_and_duplex_strand_family_sizes(run_fgumi):
    module = run_fgumi({"A.family_sizes.txt": SIMPLEX_FS, "B.family_sizes.txt": DUPLEX_FS})
    assert module.samples_parsed_by_tool["consensus_families"] == {"A", "B"}
    raw = module.saved_raw_data
    assert raw["multiqc_fgumi_simplex_family_sizes"]["A"] == {"cs": {1: 50, 3: 50}, "ss": {1: 40, 3: 60}}
    assert raw["multiqc_fgumi_duplex_family_sizes"]["B"]["ds"] == {1: 30, 2: 70}
    assert {"fgumi_simplex_family_sizes", "fgumi_duplex_family_sizes"} <= set(report.plot_by_id)


def test_duplex_ab_ba_heatmaps_and_min_strand_line(run_fgumi):
    module = run_fgumi({"B.duplex_family_sizes.txt": DUPLEX_AB_BA})
    assert module.samples_parsed_by_tool["consensus_families"] == {"B"}
    # Sizes >= 20 fold into the "20+" overflow bin; the min-strand line sums counts by min(ab, ba).
    assert module.saved_raw_data["multiqc_fgumi_duplex_min_strand"]["B"] == {0: 40, 1: 30, 22: 30}
    assert "fgumi_duplex_ab_ba_B" in report.plot_by_id
    assert "fgumi_duplex_ab_ba_both_B" in report.plot_by_id


def test_per_sample_heatmaps_are_capped_for_large_runs(run_fgumi):
    # Two heatmap sections per sample would swamp a large report; above the cap only the cross-sample line stays.
    files = {f"S{i:02d}.duplex_family_sizes.txt": DUPLEX_AB_BA for i in range(11)}
    module = run_fgumi(files)
    assert len(module.samples_parsed_by_tool["consensus_families"]) == 11
    assert "fgumi_duplex_min_strand" in report.plot_by_id
    assert not [pid for pid in report.plot_by_id if str(pid).startswith("fgumi_duplex_ab_ba")]


def test_sample_names_are_html_escaped_in_section_titles(run_fgumi):
    module = run_fgumi({"<i>x.duplex_family_sizes.txt": DUPLEX_AB_BA})
    names = [section.name for section in module.sections]
    assert any("&lt;i&gt;x" in name for name in names), names
    assert not any("<i>" in name for name in names), names


def test_heatmap_axes_are_contiguous_and_capped(run_fgumi):
    import math

    import pytest

    run_fgumi({"B.duplex_family_sizes.txt": DUPLEX_AB_BA})
    dataset = report.plot_by_id["fgumi_duplex_ab_ba_B"].datasets[0]
    assert dataset.ycats == [str(i) for i in range(1, 20)] + ["20+"]
    assert dataset.xcats == [str(i) for i in range(0, 20)] + ["20+"]
    assert dataset.rows[-1][-1] == pytest.approx(math.log10(30))


def test_both_strands_heatmap_section_kept_with_alert_when_empty(run_fgumi):
    only_ab = "ab_size\tba_size\tcount\tfraction\tfraction_gt_or_eq_size\n1\t0\t40\t0.4\t1\n2\t0\t60\t0.6\t0.6\n"
    module = run_fgumi({"B.duplex_family_sizes.txt": only_ab})
    both = [s for s in module.sections if s.anchor == "fgumi-duplex-ab-ba-both-B"]
    assert len(both) == 1 and both[0].plot_anchor is None
    assert both[0].alerts
