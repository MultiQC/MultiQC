"""Inline-fixture unit tests for the riker module.

The fixture strings here are short slices of real riker outputs (1000 Genomes
samples HG00240 / HG03953); they're enough to exercise the parsers without
embedding multi-megabyte histograms.
"""

import pytest

from multiqc import report
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.riker import MultiqcModule

ALIGNMENT_TSV = (
    "sample\tcategory\ttotal_reads\taligned_reads\tfrac_aligned\thq_aligned_reads\thq_aligned_bases"
    "\thq_aligned_q20_bases\thq_median_mismatches\tmismatch_rate\thq_mismatch_rate\tindel_rate"
    "\tmean_read_length\tsd_read_length\tmedian_read_length\tmad_read_length\tmin_read_length"
    "\tmax_read_length\tmean_aligned_read_length\taligned_reads_in_pairs\tfrac_aligned_in_pairs"
    "\treads_improperly_paired\tfrac_reads_improperly_paired\tbad_cycles\tstrand_balance"
    "\tfrac_chimeras\tfrac_softclipped_reads\tfrac_hardclipped_reads\tmean_3prime_softclipped_bases\n"
    "HG00240\tread1\t49290193\t49272959\t0.99965\t45932523\t4543946625\t4266911175\t0.00\t0.004599"
    "\t0.004363\t0.000106\t100.00\t0.00\t100.00\t0.00\t100\t100\t98.82\t49091197\t0.99631\t1213329"
    "\t0.02462\t0\t0.50414\t0.01790\t0.01184\t0.00000\t21.59\n"
    "HG00240\tread2\t49290193\t49096151\t0.99606\t45744255\t4501617843\t4174244580\t0.00\t0.006610"
    "\t0.006364\t0.000109\t100.00\t0.00\t100.00\t0.00\t100\t100\t98.25\t49091194\t0.99990\t1213326"
    "\t0.02471\t0\t0.49571\t0.01796\t0.01744\t0.00000\t21.78\n"
    "HG00240\tpair\t98580386\t98369110\t0.99786\t91676778\t9045564468\t8441155755\t0.00\t0.005604"
    "\t0.005364\t0.000108\t100.00\t0.00\t100.00\t0.00\t100\t100\t98.54\t98182391\t0.99811\t2426655"
    "\t0.02466\t0\t0.49993\t0.01793\t0.01464\t0.00000\t21.69\n"
)

BASE_DIST_TSV = (
    "sample\tread_end\tcycle\tfrac_a\tfrac_c\tfrac_g\tfrac_t\tfrac_n\n"
    "HG00240\t1\t1\t0.28942\t0.20943\t0.20676\t0.29184\t0.00255\n"
    "HG00240\t1\t2\t0.28603\t0.21469\t0.21071\t0.28856\t0.00001\n"
    "HG00240\t2\t1\t0.28800\t0.21100\t0.20800\t0.29100\t0.00200\n"
    "HG00240\t2\t2\t0.28500\t0.21300\t0.21100\t0.29000\t0.00100\n"
)

MEAN_QUAL_TSV = (
    "sample\tcycle\tmean_quality\nHG00240\t1\t26.50\nHG00240\t2\t26.34\nHG00240\t3\t27.10\nHG00240\t4\t27.45\n"
)

QUAL_DIST_TSV = (
    "sample\tquality\tcount\tfrac_bases\n"
    "HG00240\t1\t2735244\t0.00028\n"
    "HG00240\t6\t696483380\t0.07070\n"
    "HG00240\t15\t105248183\t0.01068\n"
    "HG00240\t22\t192281327\t0.01952\n"
    "HG00240\t37\t8856200000\t0.89882\n"
)

GCBIAS_DETAIL_TSV = (
    "sample\tgc\twindows\tread_starts\treported_base_quality\tempirical_base_quality"
    "\tnormalized_coverage\terror_bar_width\n"
    "HG03953\t0\t139943\t1752\t28.66\t14.42\t0.73802\t0.01763\n"
    "HG03953\t1\t102415\t3123\t29.13\t18.94\t1.79760\t0.03217\n"
    "HG03953\t2\t98765\t4500\t29.50\t20.00\t2.10000\t0.04000\n"
)

GCBIAS_SUMMARY_TSV = (
    "sample\twindow_size\ttotal_clusters\taligned_reads\tat_dropout\tgc_dropout"
    "\tgc_0_19_normcov\tgc_20_39_normcov\tgc_40_59_normcov\tgc_60_79_normcov\tgc_80_100_normcov\n"
    "HG03953\t100\t25794918\t51627320\t1.28581\t0.78299\t1.60739\t1.00695\t0.97494\t1.02991\t3.67370\n"
)

HYBCAP_METRICS_TSV = (
    "sample\tpanel_name\tgenome_size\tbait_territory\ttarget_territory\tbait_design_efficiency"
    "\ttotal_reads\tdeduped_reads\tdeduped_reads_aligned\ttotal_bases\tbases_aligned"
    "\tdeduped_bases_aligned\ton_bait_bases\tnear_bait_bases\toff_bait_bases\ton_target_bases"
    "\ton_target_bases_from_pairs\tselected_pairs\tselected_unique_pairs\tfrac_deduped_reads"
    "\tfrac_deduped_reads_aligned\tfrac_selected_bases\tfrac_off_bait\ton_bait_vs_selected"
    "\tmean_bait_coverage\tmean_target_coverage\tmedian_target_coverage\tmax_target_coverage"
    "\tmin_target_coverage\tfrac_uncovered_targets\tfold_enrichment\tfold_80_base_penalty"
    "\tfrac_exc_dupe\tfrac_exc_mapq\tfrac_exc_overlap\tfrac_exc_baseq\tfrac_exc_off_target"
    "\tfrac_target_bases_1x\tfrac_target_bases_10x\tfrac_target_bases_20x\tfrac_target_bases_30x"
    "\tfrac_target_bases_50x\tfrac_target_bases_100x\tfrac_target_bases_250x\tfrac_target_bases_500x"
    "\tfrac_target_bases_1000x\tat_dropout\tgc_dropout\tfrac_usable_bases_on_bait"
    "\tfrac_usable_bases_on_target\ths_library_size\ths_penalty_10x\ths_penalty_20x\ths_penalty_30x"
    "\ths_penalty_40x\ths_penalty_50x\ths_penalty_100x\n"
    "HG00240\toutput_1000G_Exome.v1\t3217346917\t46106635\t46106635\t1.000000\t98580386\t91792601"
    "\t91581325\t9858038600\t10643786987\t9978299743\t4925288538\t2365011060\t3353487389\t3832936583"
    "\t3830372008\t36910710\t34275479\t0.931145\t0.997698\t0.684935\t0.315065\t0.675595\t106.82"
    "\t83.13\t69.00\t3051\t0\t0.022314\t32.29\t3.61\t0.062524\t0.076124\t0.032542\t0.053970"
    "\t0.414730\t0.966910\t0.888920\t0.824546\t0.760370\t0.629690\t0.327539\t0.028106\t0.001281"
    "\t0.000048\t0.40\t19.63\t0.499622\t0.388813\t246041533\t9.70\t10.03\t10.38\t10.77\t11.21\t14.46\n"
)

ISIZE_METRICS_TSV = (
    "sample\tpair_orientation\tread_pairs\tmean_insert_size\tstandard_deviation\tmedian_insert_size"
    "\tmedian_absolute_deviation\tmode_insert_size\tmin_insert_size\tmax_insert_size\n"
    "HG00240\tFR\t44718990\t253.40\t80.08\t241.00\t50.00\t212\t2\t247959841\n"
)

ISIZE_HISTOGRAM_TSV = (
    "sample\tinsert_size\tfr_count\trf_count\ttandem_count\n"
    "HG00240\t1\t0\t3\t19\n"
    "HG00240\t2\t18\t6\t12\n"
    "HG00240\t3\t28\t9\t6\n"
    "HG00240\t100\t1000\t10\t5\n"
    "HG00240\t200\t5000\t20\t8\n"
    "HG00240\t300\t1500\t15\t4\n"
)

WGS_METRICS_TSV = (
    "sample\tgenome_territory\tmean_coverage\tsd_coverage\tmedian_coverage\tmad_coverage"
    "\tfrac_excluded_mapq\tfrac_excluded_dupe\tfrac_excluded_unpaired\tfrac_excluded_baseq"
    "\tfrac_excluded_overlap\tfrac_excluded_capped\tfrac_excluded_total\tfrac_bases_at_1x"
    "\tfrac_bases_at_5x\tfrac_bases_at_10x\tfrac_bases_at_15x\tfrac_bases_at_20x\tfrac_bases_at_25x"
    "\tfrac_bases_at_30x\tfrac_bases_at_40x\tfrac_bases_at_50x\tfrac_bases_at_60x\tfrac_bases_at_100x"
    "\tfold_80_base_penalty\tfold_90_base_penalty\tfold_95_base_penalty\n"
    "HG03953\t3043453562\t2.06\t2.10\t2.00\t1.00\t0.05880\t0.08244\t0.00091\t0.01497\t0.00666"
    "\t0.00472\t0.16850\t0.82110\t0.07184\t0.00079\t0.00031\t0.00022\t0.00017\t0.00013\t0.00009"
    "\t0.00007\t0.00006\t0.00004\t2.06\t0.00\t0.00\n"
)

WGS_COVERAGE_TSV = (
    "sample\tdepth\tbases\tfrac_bases\tbases_at_or_above\tfrac_bases_at_or_above\n"
    "HG03953\t0\t544461140\t0.17890\t3043453562\t1.00000\n"
    "HG03953\t1\t701573680\t0.23052\t2498992422\t0.82110\n"
    "HG03953\t2\t737951772\t0.24247\t1797418742\t0.59059\n"
    "HG03953\t3\t538916709\t0.17707\t1059466970\t0.34811\n"
    "HG03953\t10\t1000\t0.00001\t100\t0.00001\n"
)


@pytest.fixture
def run_riker_module(tmp_path):
    """Factory: write the given content to a temp file and run the riker module against it."""

    def _run(filename: str, content: str):
        path = tmp_path / filename
        path.write_text(content)
        report.reset()
        report.analysis_files = [path]
        report.search_files(["riker"])
        return MultiqcModule()

    return _run


def _assert_one_sample_for_tool(module, tool: str, expected: str):
    parsed = module.samples_parsed_by_tool.get(tool, set())
    assert expected in parsed, f"{tool} did not parse {expected!r}; got {parsed!r}"


def test_alignment_metrics(run_riker_module):
    module = run_riker_module("HG00240.alignment-metrics.txt", ALIGNMENT_TSV)
    _assert_one_sample_for_tool(module, "alignment", "HG00240")


def test_basic_base_distribution(run_riker_module):
    module = run_riker_module("HG00240.base-distribution-by-cycle.txt", BASE_DIST_TSV)
    # When read_end=2 is present we split into _R1 / _R2; the per-tool sample set
    # collapses both back to the underlying sample name.
    _assert_one_sample_for_tool(module, "basic", "HG00240")


def test_basic_mean_quality(run_riker_module):
    module = run_riker_module("HG00240.mean-quality-by-cycle.txt", MEAN_QUAL_TSV)
    _assert_one_sample_for_tool(module, "basic", "HG00240")


def test_basic_quality_distribution(run_riker_module):
    module = run_riker_module("HG00240.quality-score-distribution.txt", QUAL_DIST_TSV)
    _assert_one_sample_for_tool(module, "basic", "HG00240")


def test_gcbias_detail(run_riker_module):
    module = run_riker_module("HG03953.gcbias-detail.txt", GCBIAS_DETAIL_TSV)
    _assert_one_sample_for_tool(module, "gcbias", "HG03953")


def test_gcbias_summary(run_riker_module):
    module = run_riker_module("HG03953.gcbias-summary.txt", GCBIAS_SUMMARY_TSV)
    _assert_one_sample_for_tool(module, "gcbias", "HG03953")


def test_hybcap_metrics(run_riker_module):
    module = run_riker_module("HG00240.hybcap-metrics.txt", HYBCAP_METRICS_TSV)
    _assert_one_sample_for_tool(module, "hybcap", "HG00240")


def test_isize_metrics(run_riker_module):
    module = run_riker_module("HG00240.isize-metrics.txt", ISIZE_METRICS_TSV)
    _assert_one_sample_for_tool(module, "isize", "HG00240")


def test_isize_histogram(run_riker_module):
    module = run_riker_module("HG00240.isize-histogram.txt", ISIZE_HISTOGRAM_TSV)
    _assert_one_sample_for_tool(module, "isize", "HG00240")


def test_wgs_metrics(run_riker_module):
    module = run_riker_module("HG03953.wgs-metrics.txt", WGS_METRICS_TSV)
    _assert_one_sample_for_tool(module, "wgs", "HG03953")


def test_wgs_coverage(run_riker_module):
    module = run_riker_module("HG03953.wgs-coverage.txt", WGS_COVERAGE_TSV)
    _assert_one_sample_for_tool(module, "wgs", "HG03953")


def test_no_riker_files_raises(tmp_path):
    # An unrelated file should not match any riker search pattern.
    unrelated = tmp_path / "irrelevant.txt"
    unrelated.write_text("hello world\n")
    report.reset()
    report.analysis_files = [unrelated]
    report.search_files(["riker"])
    with pytest.raises(ModuleNoSamplesFound):
        MultiqcModule()
