"""Tests for the Somalier module."""

import pytest

from multiqc import config, report
from multiqc.modules.somalier import MultiqcModule

SAMPLES_TSV = (
    "#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\n"
    "FAM1\tsampleA\t0\t0\t1\t-9\n"
    "FAM1\tsampleB\t0\t0\t2\t-9\n"
)


def _pairs_tsv(concordance_col: str) -> str:
    return (
        "#sample_a\tsample_b\trelatedness\tibs0\tibs2\t"
        f"{concordance_col}\thets_a\thets_b\thets_ab\tshared_hets\t"
        "hom_alts_a\thom_alts_b\tshared_hom_alts\tn\tx_ibs0\tx_ibs2\t"
        "expected_relatedness\n"
        "sampleA\tsampleB\t0.500\t10\t100\t0.990\t50\t50\t40\t45\t20\t20\t"
        "18\t200\t0\t5\t0.5\n"
    )


@pytest.fixture
def run_somalier(tmp_path):
    """Write samples plus pairs TSVs and run the Somalier module against them."""

    def _run(pairs_text: str):
        (tmp_path / "somalier.samples.tsv").write_text(SAMPLES_TSV)
        (tmp_path / "somalier.pairs.tsv").write_text(pairs_text)
        report.reset()
        config.reset()
        report.analysis_files = [tmp_path]
        report.search_files(["somalier"])
        config.preserve_module_raw_data = True
        return MultiqcModule()

    return _run


@pytest.mark.parametrize("concordance_col", ["concordance", "hom_concordance"])
def test_relatedness_heatmap_from_pairs_tsv(run_somalier, concordance_col):
    """Newer somalier writes `concordance` instead of `hom_concordance`.

    The pairs search pattern must match both headers so the relatedness
    heatmap is still produced. Samples.tsv is present either way, so a
    missed pairs file used to yield a report with no heatmap rather than
    ModuleNoSamplesFound.
    """
    module = run_somalier(_pairs_tsv(concordance_col))

    assert "Relatedness Heatmap" in [section.name for section in module.sections]
    assert "somalier_relatedness_heatmap_plot" in report.plot_by_id

    saved = module.saved_raw_data
    assert saved is not None
    pair = saved["multiqc_somalier"]["sampleA*sampleB"]
    assert pair["relatedness"] == 0.5
