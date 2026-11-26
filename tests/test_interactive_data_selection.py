import pytest

import multiqc
from multiqc import report


@pytest.mark.parametrize("module", ["samtools", "Samtools", "SAMTOOLS"])
def test_get_module_data_by_name(data_dir, module):
    multiqc.parse_logs(
        data_dir / "modules/samtools/flagstat/small.samtools13.flagstat.log.txt",
        run_modules=["samtools"],
        no_version_check=True,
        strict=True,
    )
    data = multiqc.get_module_data(module=module)
    assert data == multiqc.get_module_data()["samtools"]
    assert next(iter(data.values()))["total_passed"] == 5414


def test_get_module_data_key_selection(data_dir):
    multiqc.parse_logs(
        data_dir / "modules/samtools/flagstat",
        data_dir / "modules/samtools/stats",
        run_modules=["samtools"],
        no_version_check=True,
        strict=True,
    )
    all_data = multiqc.get_module_data(module="samtools")
    key = "multiqc_samtools_flagstat"
    assert len(all_data) == 2
    assert multiqc.get_module_data(module="samtools", key=key) == {key: all_data[key]}
    sample = next(iter(all_data[key]))
    assert multiqc.get_module_data(module="samtools", key=key, sample=sample) == {key: all_data[key][sample]}
    assert multiqc.get_module_data(key=key) == {"samtools": {key: all_data[key]}}
    with pytest.raises(ValueError, match="Key .* is not found"):
        multiqc.get_module_data(module="samtools", key="missing")
    with pytest.raises(ValueError, match="Module .* is not found"):
        multiqc.get_module_data(module="missing")
    assert len(report.modules) == 1


def test_get_module_data_by_custom_name_and_anchor(data_dir, tmp_path):
    cfg = tmp_path / "multiqc_config.yaml"
    cfg.write_text("module_order:\n  - samtools:\n      name: Alignment QC\n      anchor: align\n")
    multiqc.parse_logs(
        data_dir / "modules/samtools/flagstat/small.samtools13.flagstat.log.txt",
        config_files=[cfg],
        run_modules=["samtools"],
        no_version_check=True,
        strict=True,
    )
    expected = multiqc.get_module_data()["align"]
    assert multiqc.get_module_data(module="Alignment QC") == expected
    assert multiqc.get_module_data(module="alignment qc") == expected
    assert multiqc.get_module_data(module="align") == expected
