import os

import multiqc
from multiqc import report


def test_parse_logs_fn_clean_exts(data_dir):
    multiqc.reset()
    multiqc.parse_logs(
        data_dir / "modules/fastp/SAMPLE.json",
        data_dir / "modules/fastp/single_end",
        extra_fn_clean_exts=["_1", "_S10_R1_001"],
    )
    assert multiqc.list_samples() == ["SRR5442949", "smalltest"]
    assert multiqc.list_modules() == ["fastp"]


def test_parse_logs_ignore_samples(data_dir):
    multiqc.reset()
    multiqc.parse_logs(
        data_dir / "modules/quast/full_metaquast_run",
        ignore_samples=["meta_contigs_2"],
    )

    assert multiqc.list_samples() == ["meta_contigs_1"]
    assert multiqc.list_modules() == ["QUAST"]


def test_write_report(tmp_path):
    multiqc.reset()
    module = multiqc.BaseMultiqcModule(name="my-module", anchor="custom_data")
    module.add_section()
    report.modules = [module]

    multiqc.write_report(force=True, output_dir=str(tmp_path))
    assert (tmp_path / "multiqc_report.html").is_file()
    assert (tmp_path / "multiqc_data").is_dir()


def test_software_versions_section(data_dir, tmp_path, capsys):
    multiqc.reset()

    multiqc.parse_logs(data_dir / "modules/fastp")
    multiqc.parse_logs(data_dir / "modules/bcftools")
    multiqc.write_report(filename="stdout")  # triggers adding software_versions module
    assert multiqc.list_modules() == ["fastp", "Bcftools", "Software Versions"]


def test_write_report_multiple_times(data_dir, tmp_path):
    multiqc.reset()
    multiqc.parse_logs(data_dir / "modules/fastp")
    multiqc.write_report(output_dir=str(tmp_path))
    assert multiqc.list_modules() == ["fastp", "Software Versions"]

    multiqc.parse_logs(data_dir / "modules/bcftools")
    multiqc.write_report(output_dir=str(tmp_path))
    assert multiqc.list_modules() == ["fastp", "Bcftools", "Software Versions"]


def test_run_twice(data_dir, tmp_path):
    from multiqc import multiqc
    from multiqc.core.update_config import ClConfig

    data_dir = data_dir / "custom_content/with_config/run_concordance"
    os.chdir(tmp_path)
    multiqc.run(  # pylint: disable=no-member
        data_dir / "run_concordance.txt",
        cfg=ClConfig(
            strict=True,
            force=True,
            config_files=[data_dir / "multiqc_config.yaml"],
        ),
    )
    assert (tmp_path / "multiqc_report.html").is_file()
    assert (tmp_path / "multiqc_data").is_dir()

    multiqc.run(  # pylint: disable=no-member
        data_dir / "run_concordance.txt",
        cfg=ClConfig(
            strict=True,
            force=True,
            config_files=[data_dir / "multiqc_config.yaml"],
        ),
    )
    assert (tmp_path / "multiqc_report.html").is_file()
    assert (tmp_path / "multiqc_data").is_dir()
