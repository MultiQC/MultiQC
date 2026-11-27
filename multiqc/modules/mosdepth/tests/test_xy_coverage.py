import pytest

from multiqc import config, report
from multiqc.modules.mosdepth import MultiqcModule
from multiqc.plots.bargraph import BarPlot
from multiqc.plots.linegraph import LinePlot
from multiqc.types import Anchor

# Three rows per contig: the cumulative fractions mosdepth writes. Mean coverage
# is the sum over depth >= 1 (excluding depth 0, which is always 1.00), so the X
# contig comes to 1.0 and the Y contig to 0.3 whatever the contigs are called.
DIST_TEMPLATE = """\
chr1\t2\t0.90
chr1\t1\t0.80
chr1\t0\t1.00
{xchr}\t2\t0.60
{xchr}\t1\t0.40
{xchr}\t0\t1.00
{ychr}\t2\t0.20
{ychr}\t1\t0.10
{ychr}\t0\t1.00
total\t2\t0.70
total\t1\t0.50
total\t0\t1.00
"""


@pytest.mark.parametrize(
    "x_contig,y_contig,mosdepth_config,expected_names",
    [
        ("myXchr", "myYchr", {"xchr": "myXchr", "ychr": "myYchr"}, ["myXchr", "myYchr"]),
        ("myXchr", "chrY", {"xchr": "myXchr"}, ["myXchr", "Chromosome Y"]),
        ("chrX", "myYchr", {"ychr": "myYchr"}, ["Chromosome X", "myYchr"]),
        ("chrX", "chrY", {}, ["Chromosome X", "Chromosome Y"]),
    ],
    ids=["both_keys", "xchr_only", "ychr_only", "no_keys"],
)
def test_xy_coverage_bar_names(monkeypatch, tmp_path, x_contig, y_contig, mosdepth_config, expected_names):
    """Each bar of the XY coverage plot is named from its own config key.

    The X and Y coverage are collected separately, so the names have to line up
    with the values: a bar named for the X chromosome must carry the X coverage.
    """
    dist_file = tmp_path / "sampleA.mosdepth.global.dist.txt"
    dist_file.write_text(DIST_TEMPLATE.format(xchr=x_contig, ychr=y_contig))

    # mosdepth_config is not a config default, so config.reset() would not clear
    # it and it would leak into every later test in the process.
    monkeypatch.setattr(config, "mosdepth_config", dict(mosdepth_config), raising=False)

    report.reset()
    report.analysis_files = [dist_file]
    report.search_files(["mosdepth"])
    MultiqcModule()

    plot = report.plot_by_id[Anchor("mosdepth-xy-coverage-plot")]
    assert isinstance(plot, BarPlot)
    cats = plot.datasets[0].cats
    assert [cat.name for cat in cats] == expected_names
    assert [list(cat.data) for cat in cats] == [[pytest.approx(1.0)], [pytest.approx(0.3)]]


def test_xy_coverage_matches_per_contig_with_zero_y(monkeypatch, tmp_path):
    """XY coverage excludes the depth-0 row so it matches per-contig average coverage, including 0X chrY."""
    female_file = tmp_path / "female.mosdepth.global.dist.txt"
    female_file.write_text(
        """\
chr1\t2\t0.90
chr1\t1\t0.95
chr1\t0\t1.00
chrX\t2\t0.60
chrX\t1\t0.80
chrX\t0\t1.00
chrY\t0\t1.00
total\t2\t0.70
total\t1\t0.85
total\t0\t1.00
"""
    )
    male_file = tmp_path / "male.mosdepth.global.dist.txt"
    male_file.write_text(
        """\
chr1\t2\t0.90
chr1\t1\t0.95
chr1\t0\t1.00
chrX\t2\t0.30
chrX\t1\t0.40
chrX\t0\t1.00
chrY\t2\t0.20
chrY\t1\t0.30
chrY\t0\t1.00
total\t2\t0.60
total\t1\t0.70
total\t0\t1.00
"""
    )
    monkeypatch.setattr(config, "mosdepth_config", {}, raising=False)

    report.reset()
    report.analysis_files = [female_file, male_file]
    report.search_files(["mosdepth"])
    MultiqcModule()

    xy_plot = report.plot_by_id[Anchor("mosdepth-xy-coverage-plot")]
    assert isinstance(xy_plot, BarPlot)
    ds = xy_plot.datasets[0]
    xy_by_sample = {sample: {cat.name: cat.data[i] for cat in ds.cats} for i, sample in enumerate(ds.samples)}

    perchrom_plot = report.plot_by_id[Anchor("mosdepth-coverage-per-contig-multi")]
    assert isinstance(perchrom_plot, LinePlot)
    perchrom_by_sample = {line.name: dict(line.pairs) for line in perchrom_plot.datasets[0].lines}

    assert (
        xy_by_sample["female"]["Chromosome X"]
        == pytest.approx(perchrom_by_sample["female"]["chrX"])
        == pytest.approx(1.4)
    )
    assert (
        xy_by_sample["female"]["Chromosome Y"]
        == pytest.approx(perchrom_by_sample["female"]["chrY"])
        == pytest.approx(0.0)
    )
    assert (
        xy_by_sample["male"]["Chromosome X"] == pytest.approx(perchrom_by_sample["male"]["chrX"]) == pytest.approx(0.7)
    )
    assert (
        xy_by_sample["male"]["Chromosome Y"] == pytest.approx(perchrom_by_sample["male"]["chrY"]) == pytest.approx(0.5)
    )


def test_perchrom_fraction_cutoff_excludes_zero_coverage_contig(monkeypatch, tmp_path):
    """perchrom_fraction_cutoff filters on mean coverage without adding 1.0 from depth-0 rows."""
    dist_file = tmp_path / "sampleA.mosdepth.global.dist.txt"
    dist_file.write_text(
        """\
chr1\t2\t0.90
chr1\t1\t0.90
chr1\t0\t1.00
chr2\t0\t1.00
total\t2\t0.45
total\t1\t0.45
total\t0\t1.00
"""
    )
    monkeypatch.setattr(config, "mosdepth_config", {"perchrom_fraction_cutoff": 0.1}, raising=False)

    report.reset()
    report.analysis_files = [dist_file]
    report.search_files(["mosdepth"])
    MultiqcModule()

    # Only chr1 passes the cutoff, so the single-contig bar plot is emitted
    single_plot = report.plot_by_id[Anchor("mosdepth-coverage-per-contig-single")]
    assert isinstance(single_plot, BarPlot)
    assert [cat.name for cat in single_plot.datasets[0].cats] == ["chr1"]
    assert single_plot.datasets[0].cats[0].data[0] == pytest.approx(1.8)
