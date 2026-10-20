import gzip

import pytest

from multiqc import config, report, reset
from multiqc.base_module import ModuleNoSamplesFound


@pytest.fixture(autouse=True)
def _reset_report():
    """Isolate from file-search state left behind by other tests running in the same process."""
    reset()


def test_thresholds_bed_parse_error_is_skipped_not_crashed(tmp_path, caplog):
    """A malformed thresholds.bed.gz (e.g. a same-named file from another tool) should be skipped
    with a warning, not crash the module - *.thresholds.bed.gz is a generic enough name that a
    coincidental match from an unrelated tool is expected to happen."""
    thresholds_file = tmp_path / "sample.thresholds.bed.gz"
    thresholds_file.write_bytes(gzip.compress(b"chrom\tstart\tend\tname\t20X\n1\t100\t200\tGENE1\t50\n"))

    report.analysis_files = [tmp_path]
    report.search_files(["mosdepth"])

    from multiqc.modules.mosdepth.mosdepth import MultiqcModule

    with pytest.raises(ModuleNoSamplesFound):
        MultiqcModule()

    assert "sample.thresholds.bed.gz" in caplog.text


@pytest.mark.parametrize(
    "content",
    [
        b"1\t100\t200\tGENE1\t45.5\n",  # plain text, not gzip
        gzip.compress(b"1\t100\t200\tGENE1\t45.5\n" * 100)[:-20],  # truncated
    ],
    ids=["not-gzip", "truncated"],
)
def test_unreadable_regions_bed_is_skipped_not_crashed(tmp_path, caplog, content):
    (tmp_path / "sample.regions.bed.gz").write_bytes(content)

    report.analysis_files = [tmp_path]
    report.search_files(["mosdepth"])

    from multiqc.modules.mosdepth.mosdepth import MultiqcModule

    with pytest.raises(ModuleNoSamplesFound):
        MultiqcModule()

    assert "sample.regions.bed.gz" in caplog.text


def test_regions_bed_over_max_table_rows_is_skipped(tmp_path, caplog, monkeypatch):
    monkeypatch.setattr(config, "max_table_rows", 2)
    rows = b"".join(b"1\t%d\t%d\t1.0\n" % (i * 10, i * 10 + 10) for i in range(3))
    (tmp_path / "sample.regions.bed.gz").write_bytes(gzip.compress(rows))

    report.analysis_files = [tmp_path]
    report.search_files(["mosdepth"])

    from multiqc.modules.mosdepth.mosdepth import MultiqcModule

    with pytest.raises(ModuleNoSamplesFound):
        MultiqcModule()

    assert "More than 2 regions" in caplog.text


def test_regions_bed_without_thresholds_gets_mean_coverage_only_section(tmp_path, monkeypatch):
    """regions.bed.gz alone (--by without --thresholds) should still add the per-region section,
    with mean coverage but no threshold percentage columns."""
    regions_file = tmp_path / "sample.regions.bed.gz"
    regions_file.write_bytes(gzip.compress(b"1\t100\t200\tGENE1\t45.5\n"))

    monkeypatch.setattr(config, "preserve_module_raw_data", True)
    report.analysis_files = [tmp_path]
    report.search_files(["mosdepth"])

    from multiqc.modules.mosdepth.mosdepth import MultiqcModule

    module = MultiqcModule()

    assert "mosdepth-per-region-coverage" in {section.anchor for section in module.sections}
    assert module.saved_raw_data is not None
    row = module.saved_raw_data["mosdepth_per_region_coverage"]["sample | GENE1"]
    assert row == {"coordinates": "1:100-200", "mean_coverage": 45.5}
