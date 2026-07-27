import gzip

import pytest

from multiqc import report, reset


@pytest.fixture(autouse=True)
def _reset_report():
    """Isolate from file-search state left behind by other tests running in the same process."""
    reset()


def test_thresholds_bed_parse_error_includes_filename(tmp_path):
    """A malformed thresholds.bed.gz should raise with the offending filename, not a bare traceback."""
    thresholds_file = tmp_path / "sample.thresholds.bed.gz"
    thresholds_file.write_bytes(gzip.compress(b"chrom\tstart\tend\tname\t20X\n1\t100\t200\tGENE1\t50\n"))

    report.analysis_files = [tmp_path]
    report.search_files(["mosdepth"])

    from multiqc.modules.mosdepth.mosdepth import MultiqcModule

    with pytest.raises(ValueError, match="sample.thresholds.bed.gz"):
        MultiqcModule()


def test_regions_bed_without_thresholds_skips_per_region_section(tmp_path):
    """regions.bed.gz alone (no --thresholds run) should be picked up, but not add the per-region section."""
    regions_file = tmp_path / "sample.regions.bed.gz"
    regions_file.write_bytes(gzip.compress(b"1\t100\t200\tGENE1\t45.5\n"))

    report.analysis_files = [tmp_path]
    report.search_files(["mosdepth"])

    from multiqc.modules.mosdepth.mosdepth import MultiqcModule

    module = MultiqcModule()

    assert "mosdepth-per-region-coverage" not in {section.anchor for section in module.sections}
