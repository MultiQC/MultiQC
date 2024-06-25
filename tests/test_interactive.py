import tempfile
from pathlib import Path

import multiqc
from multiqc import report, reset


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


def test_write_report():
    multiqc.reset()
    module = multiqc.BaseMultiqcModule(
        name="my-module",
        anchor="custom_data",
    )
    module.add_section()
    report.modules = [module]

    with tempfile.TemporaryDirectory() as tmp_dir_name:
        multiqc.write_report(force=True, output_dir=tmp_dir_name)
        assert (Path(tmp_dir_name) / "multiqc_report.html").is_file()
        assert (Path(tmp_dir_name) / "multiqc_data").is_dir()


def test_run_twice(data_dir):
    from multiqc import multiqc
    from multiqc.core.update_config import ClConfig

    with tempfile.TemporaryDirectory() as tmp_dir_name:
        multiqc.run(  # pylint: disable=no-member
            data_dir / "custom_content/with_config/run_concordance/run_concordance.txt",
            cfg=ClConfig(
                strict=True,
                force=True,
                config_files=[data_dir / "custom_content/with_config/run_concordance/multiqc_config.yaml"],
                output_dir=tmp_dir_name,
            ),
        )
        path = Path(tmp_dir_name) / "multiqc_report.html"
        assert path.is_file()
        path = Path(tmp_dir_name) / "multiqc_data"
        assert path.is_dir()

    with tempfile.TemporaryDirectory() as tmp_dir_name:
        multiqc.run(  # pylint: disable=no-member
            data_dir / "custom_content/with_config/run_concordance/run_concordance.txt",
            cfg=ClConfig(
                strict=True,
                force=True,
                config_files=[data_dir / "custom_content/with_config/run_concordance/multiqc_config.yaml"],
                output_dir=tmp_dir_name,
            ),
        )
        path = Path(tmp_dir_name) / "multiqc_report.html"
        assert path.is_file()
        path = Path(tmp_dir_name) / "multiqc_data"
        assert path.is_dir()
