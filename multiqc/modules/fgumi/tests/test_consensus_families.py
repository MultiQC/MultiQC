import pytest

from multiqc import config, report

from .conftest import line_points

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
    assert raw["multiqc_fgumi_simplex_family_sizes"]["A"] == {"cs_1": 50, "cs_3": 50, "ss_1": 40, "ss_3": 60}
    duplex = raw["multiqc_fgumi_duplex_family_sizes"]["B"]
    assert duplex["ds_1"] == 30
    assert duplex["ds_2"] == 70
    assert {"fgumi_simplex_family_sizes", "fgumi_duplex_family_sizes"} <= set(report.plot_by_id)


def test_duplex_ab_ba_heatmaps_and_min_strand_line(run_fgumi):
    module = run_fgumi({"B.duplex_family_sizes.txt": DUPLEX_AB_BA})
    assert module.samples_parsed_by_tool["consensus_families"] == {"B"}
    # Sizes >= 20 fold into the "20+" overflow bin; the min-strand line sums counts by min(ab, ba).
    assert module.saved_raw_data["multiqc_fgumi_duplex_min_strand"]["B"] == {0: 40, 1: 30, 22: 30}
    assert "fgumi_duplex_ab_ba_B" in report.plot_by_id
    assert "fgumi_duplex_ab_ba_both_B" in report.plot_by_id


@pytest.mark.parametrize(
    "n_samples, cap, expect_heatmaps",
    [
        pytest.param(6, None, False, id="above-default-cap"),
        pytest.param(5, None, True, id="at-default-cap"),
        pytest.param(2, 1, False, id="configured-cap"),
    ],
)
def test_per_sample_heatmaps_are_capped_for_large_runs(run_fgumi, n_samples, cap, expect_heatmaps):
    # Two heatmap sections per sample would swamp a large report; above the cap only the cross-sample line stays.
    if cap is not None:
        config.fgumi_config = {"max_duplex_heatmap_samples": cap}
    files = {f"S{i:02d}.duplex_family_sizes.txt": DUPLEX_AB_BA for i in range(n_samples)}
    module = run_fgumi(files)
    assert len(module.samples_parsed_by_tool["consensus_families"]) == n_samples
    assert "fgumi_duplex_min_strand" in report.plot_by_id
    heatmaps = [pid for pid in report.plot_by_id if str(pid).startswith("fgumi_duplex_ab_ba")]
    assert bool(heatmaps) == expect_heatmaps


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
    # As fgbio draws it: AB size on x, BA size on y.
    assert dataset.xcats == [str(i) for i in range(1, 20)] + ["20+"]
    assert dataset.ycats == [str(i) for i in range(0, 20)] + ["20+"]
    assert dataset.rows[-1][-1] == pytest.approx(math.log10(30))


def test_both_strands_heatmap_section_kept_with_alert_when_empty(run_fgumi):
    only_ab = "ab_size\tba_size\tcount\tfraction\tfraction_gt_or_eq_size\n1\t0\t40\t0.4\t1\n2\t0\t60\t0.6\t0.6\n"
    module = run_fgumi({"B.duplex_family_sizes.txt": only_ab})
    both = [s for s in module.sections if s.anchor == "fgumi-duplex-ab-ba-both-B"]
    assert len(both) == 1 and both[0].plot_anchor is None
    assert both[0].alerts


def test_non_finite_cumulative_points_are_left_out(run_fgumi):
    run_fgumi({"A.family_sizes.txt": SIMPLEX_FS, "B.family_sizes.txt": DUPLEX_FS})
    # Tabs: counts per strand, then cumulative per strand. NaN / Infinity cells are absent, never plotted as 0.
    assert line_points("fgumi_simplex_family_sizes", 3, "A") == {1: 100.0}  # ss cumulative; size 3 is NaN
    assert line_points("fgumi_duplex_family_sizes", 5, "B") == {1: 100.0}  # ds cumulative; size 2 is Infinity
    assert line_points("fgumi_simplex_family_sizes", 2, "A") == {1: 100.0, 3: 50.0}  # cs cumulative


def test_header_only_duplex_family_sizes_gets_one_section(run_fgumi):
    module = run_fgumi({"B.duplex_family_sizes.txt": DUPLEX_AB_BA.splitlines()[0] + "\n"})
    per_sample = [s for s in module.sections if s.anchor.startswith("fgumi-duplex-ab-ba")]
    assert [s.anchor for s in per_sample] == ["fgumi-duplex-ab-ba-B"]
    assert "No duplex families were recorded" in per_sample[0].alerts[0].message


def test_both_strands_heatmap_is_skipped_when_identical(run_fgumi):
    # Every family has BA reads, so the both-strands heatmap would repeat the all-families one.
    all_have_ba = "ab_size\tba_size\tcount\tfraction\tfraction_gt_or_eq_size\n2\t1\t30\t0.5\t1\n3\t2\t30\t0.5\t0.5\n"
    run_fgumi({"B.duplex_family_sizes.txt": all_have_ba})
    assert "fgumi_duplex_ab_ba_B" in report.plot_by_id
    assert "fgumi_duplex_ab_ba_both_B" not in report.plot_by_id


def test_duplex_metrics_files_are_all_listed_as_sources(run_fgumi):
    # duplex-metrics writes both <prefix>.family_sizes.txt and <prefix>.duplex_family_sizes.txt for one sample.
    run_fgumi({"B.family_sizes.txt": DUPLEX_FS, "B.duplex_family_sizes.txt": DUPLEX_AB_BA})
    sections = [section for section, by_sample in report.data_sources["fgumi"].items() if "B" in by_sample]
    assert len(sections) == 2
