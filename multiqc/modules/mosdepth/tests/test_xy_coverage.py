import pytest

from multiqc import config, report
from multiqc.modules.mosdepth import MultiqcModule
from multiqc.plots.bargraph import BarPlot
from multiqc.types import Anchor

# Three rows per contig: the cumulative fractions mosdepth writes. The module
# sums them, so the X contig comes to 2.0 and the Y contig to 1.3 whatever the
# contigs are called.
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
    assert [list(cat.data) for cat in cats] == [[2.0], [1.3]]
