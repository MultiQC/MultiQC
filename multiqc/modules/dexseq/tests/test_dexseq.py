import pytest

from multiqc import config, report
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.dexseq.counts import parse_counts_file

COUNTS = (
    "ENSG00000004961:001\t0\n"
    "ENSG00000004961:002\t5\n"
    "ENSG00000004961:003\t7\n"
    "ENSG00000207165+ENSG00000147403:001\t12\n"
    "_ambiguous\t56\n"
    "_ambiguous_readpair_position\t3\n"
    "_empty\t76\n"
    "_lowaqual\t2\n"
    "_notaligned\t1\n"
)


@pytest.fixture
def run_dexseq(tmp_path):
    """Write ``{filename: content}`` to a temp dir and run the dexseq module on it."""
    from multiqc.modules.dexseq import MultiqcModule

    config.reset()

    def _run(files):
        paths = []
        for name, content in files.items():
            path = tmp_path / name
            path.write_text(content)
            paths.append(path)
        report.reset()
        config.preserve_module_raw_data = True
        report.analysis_files = paths
        report.search_files(["dexseq"])
        return MultiqcModule()

    yield _run
    config.reset()


def test_counts_full():
    data = parse_counts_file(COUNTS.splitlines())
    assert data == {
        "ambiguous": 56,
        "ambiguous_readpair_position": 3,
        "empty": 76,
        "lowaqual": 2,
        "notaligned": 1,
        "bins_total": 4,
        "bins_detected": 3,
        "bins_not_detected": 1,
        "bin_counts": 24,
        "unassigned": 138,
        "bins_detected_pct": 75.0,
    }


def test_counts_without_readpair_counter():
    lines = [line for line in COUNTS.splitlines() if not line.startswith("_ambiguous_readpair")]
    data = parse_counts_file(lines)
    assert data is not None
    assert "ambiguous_readpair_position" not in data
    assert data["unassigned"] == 135


def test_counts_quoted_ids():
    lines = ['"ENSG1":"001"\t4', '"ENSG1":"002"\t0'] + [line for line in COUNTS.splitlines() if line.startswith("_")]
    data = parse_counts_file(lines)
    assert data is not None
    assert data["bins_total"] == 2
    assert data["bins_detected"] == 1


def test_counts_not_dexseq():
    assert parse_counts_file(["GENE1:001\t4", "GENE2:001\t5"]) is None


def test_counts_missing_counter():
    lines = [line for line in COUNTS.splitlines() if not line.startswith("_lowaqual")]
    with pytest.raises(ValueError, match="_lowaqual"):
        parse_counts_file(lines)


def test_counts_malformed_line():
    with pytest.raises(ValueError, match="line 2"):
        parse_counts_file(["GENE1:001\t4", "GENE1:002\t4\textra"])


def test_counts_no_bins():
    lines = [line for line in COUNTS.splitlines() if line.startswith("_")]
    data = parse_counts_file(lines)
    assert data is not None
    assert data["bins_total"] == 0
    assert "bins_detected_pct" not in data


def test_module_run(run_dexseq):
    module = run_dexseq(
        {
            "sample_a.txt": COUNTS,
            "sample_b.txt": COUNTS.replace("_empty\t76", "_empty\t10"),
            "sample_c.txt": COUNTS.replace("_lowaqual\t2\n", ""),
        }
    )
    raw = module.saved_raw_data
    assert set(raw["multiqc_dexseq_counts"]) == {"sample_a", "sample_b"}
    assert raw["multiqc_dexseq_counts"]["sample_b"]["empty"] == 10
    assert [section.anchor for section in module.sections] == ["dexseq_exon_bin_detection"]


def test_module_ignores_htseq(run_dexseq):
    htseq = "ENSG00000000003\t10\nENSG00000000005\t0\n__no_feature\t5\n__ambiguous\t1\n"
    with pytest.raises(ModuleNoSamplesFound):
        run_dexseq({"htseq.txt": htseq})
