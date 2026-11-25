import gzip
import io

import pytest

from multiqc import config, report
from multiqc.modules.leafcutter import MultiqcModule
from multiqc.modules.leafcutter.clustering import parse_perind_numers
from multiqc.modules.leafcutter.differential_splicing import parse_cluster_significance
from multiqc.types import ColumnKey, SampleGroup, SectionKey

# leafcutter-ds and leafcutter2 header order, with the optional annotations column
SIGNIFICANCE_PY = """cluster\tstatus\tloglr\tdf\tp\tp.adjust\tgenes\tannotations
chr1:clu_1_+\tSuccess\t12.5\t2.0\t1e-06\t0.0001\tGENE1\tPR
chr1:clu_2_+\tSuccess\t0.3\t1.0\t0.4\t0.6\tGENE2\tPR
chr1:clu_3_-\tSuccess\tNA\tNA\tNA\tNA\tGENE3\tNE
chr1:clu_4_-\tToo many introns in cluster\tNA\tNA\tNA\tNA\tGENE4\tPR
chr1:clu_5_-\t<=1 junction in cluster\tNA\tNA\tNA\tNA\tGENE5\tPR
chr1:clu_6_-\t<=1 sample with coverage>0\tNA\tNA\tNA\tNA\tNA\tPR
chr1:clu_7_-\t<=1 sample with coverage>min_coverage\tNA\tNA\tNA\tNA\tGENE7\tPR
chr1:clu_8_-\t<2 introns used in >=min_samples_per_intron samples\tNA\tNA\tNA\tNA\tGENE8\tPR
chr1:clu_9_-\tNot enough valid samples\tNA\tNA\tNA\tNA\tGENE9\tPR
"""

# leafcutter_ds.R without --exon_file: cluster_results_table() adds `cluster` after `p`
SIGNIFICANCE_R_NO_EXONS = """status\tloglr\tdf\tp\tcluster\tp.adjust
Success\t20.1\t3\t1e-09\tchr2:clu_1_NA\t1e-07
timeout\tNA\tNA\tNA\tchr2:clu_2_NA\tNA
Error in solve.default(H): Lapack routine dgesv: system is exactly singular\tNA\tNA\tNA\tchr2:clu_3_NA\tNA
"""

PERIND_NUMERS = """s1 s2 s3
chr1:100:200:clu_1_+ 10 0 3
chr1:100:300:clu_1_+ 5 0 0
chr1:400:500:clu_2_- 0 0 7
chr1:400:600:clu_2_- 1 0 0
"""


def test_cluster_significance_all_statuses():
    counts = parse_cluster_significance(io.StringIO(SIGNIFICANCE_PY), 0.05)
    assert counts == {
        "significant": 1,
        "not_significant": 1,
        "no_pvalue": 1,
        "timeout": 0,
        "error": 0,
        "too_many_introns": 1,
        "one_junction": 1,
        "one_sample_covered": 1,
        "one_sample_min_coverage": 1,
        "too_few_introns_used": 1,
        "not_enough_samples": 1,
        "clusters": 9,
        "tested": 3,
    }


def test_cluster_significance_r_column_order_timeout_and_error():
    counts = parse_cluster_significance(io.StringIO(SIGNIFICANCE_R_NO_EXONS), 0.05)
    assert counts["significant"] == 1
    assert counts["timeout"] == 1
    assert counts["error"] == 1
    assert counts["clusters"] == 3


def test_cluster_significance_threshold():
    counts = parse_cluster_significance(io.StringIO(SIGNIFICANCE_PY), 0.00001)
    assert counts["significant"] == 0
    assert counts["not_significant"] == 2


def test_cluster_significance_missing_column():
    with pytest.raises(ValueError, match="p.adjust"):
        parse_cluster_significance(io.StringIO("cluster\tstatus\tp\nchr1:clu_1\tSuccess\t0.1\n"), 0.05)


def test_perind_numers():
    parsed = parse_perind_numers(io.StringIO(PERIND_NUMERS))
    assert parsed["s1"] == {
        "junction_reads": 16,
        "introns": 4,
        "introns_detected": 3,
        "introns_not_detected": 1,
        "introns_detected_pct": pytest.approx(75.0),
        "clusters": 2,
        "clusters_detected": 2,
    }
    assert parsed["s2"]["junction_reads"] == 0
    assert parsed["s2"]["clusters_detected"] == 0
    assert parsed["s3"]["clusters_detected"] == 2


def test_perind_numers_without_strand_suffix():
    parsed = parse_perind_numers(io.StringIO("s1\nchr1:100:200:clu_1 4\n"))
    assert parsed["s1"]["junction_reads"] == 4


def test_perind_numers_header_only():
    parsed = parse_perind_numers(io.StringIO("s1 s2\n"))
    assert parsed["s1"]["introns"] == 0
    assert "introns_detected_pct" not in parsed["s1"]


@pytest.mark.parametrize(
    "content",
    [
        "s1 s2\nchr1:100:200:clu_1_+ 4\n",
        "s1\nnot_an_intron 4\n",
        "s1\nchr1:100:200:clu_1_+ four\n",
        "",
    ],
)
def test_perind_numers_malformed(content):
    with pytest.raises(ValueError):
        parse_perind_numers(io.StringIO(content))


@pytest.fixture
def run_module(tmp_path):
    def _run(files: dict[str, str | bytes], module_config: dict | None = None):
        for name, content in files.items():
            if isinstance(content, bytes):
                (tmp_path / name).write_bytes(content)
            elif name.endswith(".gz"):
                with gzip.open(tmp_path / name, "wt") as fh:
                    fh.write(content)
            else:
                (tmp_path / name).write_text(content)
        report.reset()
        config.reset()
        config.preserve_module_raw_data = True
        if module_config is not None:
            config.update({"leafcutter": module_config})
        report.analysis_files = [tmp_path]
        report.search_files(["leafcutter"])
        return MultiqcModule()

    return _run


def test_module_all_files(run_module):
    module = run_module(
        {
            "study_perind_numers.counts.gz": PERIND_NUMERS,
            "treated_vs_control_cluster_significance.txt": SIGNIFICANCE_PY,
        }
    )
    assert [s.anchor for s in module.sections] == [
        "leafcutter_introns_detected",
        "leafcutter_cluster_test_outcomes",
    ]
    general_stats = report.general_stats_data[SectionKey("leafcutter")][SampleGroup("s1")][0].data
    assert general_stats[ColumnKey("introns_detected_pct")] == pytest.approx(75.0)
    assert general_stats[ColumnKey("junction_reads")] == 16
    # Comparisons are not samples
    assert SampleGroup("treated_vs_control") not in report.general_stats_data[SectionKey("leafcutter")]


def test_module_r_format_without_exon_file(run_module):
    module = run_module({"ctrl_vs_kd_cluster_significance.txt": SIGNIFICANCE_R_NO_EXONS})
    assert [s.anchor for s in module.sections] == ["leafcutter_cluster_test_outcomes"]


@pytest.mark.parametrize("module_config, n_significant", [(None, 1), ({"padj_cutoff": 0.00001}, 0)])
def test_module_padj_cutoff_config(run_module, module_config, n_significant):
    run_module({"ctrl_vs_kd_cluster_significance.txt": SIGNIFICANCE_PY}, module_config)
    data = report.saved_raw_data["multiqc_leafcutter_cluster_significance"]
    assert data["ctrl_vs_kd"]["significant"] == n_significant


def test_module_skips_truncated_gzip(run_module):
    data = gzip.compress(PERIND_NUMERS.encode())
    run_module(
        {
            "study_perind_numers.counts.gz": data[: len(data) // 2],
            "ctrl_vs_kd_cluster_significance.txt": SIGNIFICANCE_PY,
            "bad_cluster_significance.txt": SIGNIFICANCE_PY.replace("\tPR\n", "\n", 1),
        }
    )
    assert "multiqc_leafcutter_clustering" not in report.saved_raw_data
    assert list(report.saved_raw_data["multiqc_leafcutter_cluster_significance"]) == ["ctrl_vs_kd"]
