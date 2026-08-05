import gzip

import pytest

from multiqc import report, reset
from multiqc.base_module import ModuleNoSamplesFound


@pytest.fixture(autouse=True)
def _reset_report():
    """Isolate from file-search state left behind by other tests running in the same process."""
    reset()


def test_thresholds_bed_parse_error_is_skipped_not_crashed(tmp_path, caplog):
    """A malformed thresholds.bed.gz (e.g. a same-named file from another tool) should be skipped
    with a debug log, not crash the module - *.thresholds.bed.gz is a generic enough name that a
    coincidental match from an unrelated tool is expected to happen."""
    thresholds_file = tmp_path / "sample.thresholds.bed.gz"
    thresholds_file.write_bytes(gzip.compress(b"chrom\tstart\tend\tname\t20X\n1\t100\t200\tGENE1\t50\n"))

    report.analysis_files = [tmp_path]
    report.search_files(["mosdepth"])

    from multiqc.modules.mosdepth.mosdepth import MultiqcModule

    with caplog.at_level("DEBUG", logger="multiqc.modules.mosdepth.mosdepth"):
        with pytest.raises(ModuleNoSamplesFound):
            MultiqcModule()

    assert "sample.thresholds.bed.gz" in caplog.text


def test_regions_bed_without_thresholds_skips_per_region_section(tmp_path):
    """regions.bed.gz alone (no --thresholds run) should be picked up, but not add the per-region section."""
    regions_file = tmp_path / "sample.regions.bed.gz"
    regions_file.write_bytes(gzip.compress(b"1\t100\t200\tGENE1\t45.5\n"))

    report.analysis_files = [tmp_path]
    report.search_files(["mosdepth"])

    from multiqc.modules.mosdepth.mosdepth import MultiqcModule

    module = MultiqcModule()

    assert "mosdepth-per-region-coverage" not in {section.anchor for section in module.sections}
