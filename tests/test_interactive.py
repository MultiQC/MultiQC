import os

import pytest

import multiqc
from multiqc import report
from multiqc.core.update_config import ClConfig
from multiqc.plots import table


def test_multiqc_run(data_dir, tmp_path):
    """
    Verify HTML and data directory with default names are written to current dir
    """
    os.chdir(tmp_path)

    files_before = set(os.listdir(os.getcwd()))

    multiqc.run(
        data_dir / "modules" / "fastp" / "SAMPLE.json",
        cfg=ClConfig(run_modules=["fastp"]),
    )

    files_after = set(os.listdir(os.getcwd()))
    assert files_after - files_before == {"multiqc_report.html", "multiqc_data"}


@pytest.mark.parametrize("clean_up", [True, False])
def test_multiqc_run_no_analysis_found(monkeypatch, tmp_path, clean_up):
    """
    Verify that an error is raised when a module is not found
    """
    import tempfile

    report_tmp_dir = tmp_path / "report_tmp"
    report_tmp_dir.mkdir()
    monkeypatch.setattr(tempfile, "mkdtemp", lambda: report_tmp_dir)

    nonexistent_file = tmp_path / "nonexistent"
    result = multiqc.run(nonexistent_file, clean_up=clean_up)

    assert result.sys_exit_code == 1
    assert result.message == "No analysis results found"
    if clean_up:
        assert not report_tmp_dir.exists()
    else:
        assert report_tmp_dir.exists()


def test_parse_logs_fn_clean_exts(data_dir):
    multiqc.parse_logs(
        data_dir / "modules/fastp/SAMPLE.json",
        data_dir / "modules/fastp/single_end",
        extra_fn_clean_exts=["_1", "_S10_R1_001"],
    )
    assert multiqc.list_samples() == ["SRR5442949", "smalltest"]
    assert multiqc.list_modules() == ["fastp"]


def test_parse_logs_ignore_samples(data_dir):
    multiqc.parse_logs(
        data_dir / "modules/quast/full_metaquast_run",
        ignore_samples=["meta_contigs_2"],
    )

    assert multiqc.list_samples() == ["meta_contigs_1"]
    assert multiqc.list_modules() == ["quast"]


def test_custom_module(tmp_path):
    module = multiqc.BaseMultiqcModule(name="my-module", anchor="custom_data")
    module.add_section(
        name="Custom Section",
        description="Custom description",
        helptext="Custom help",
        plot=table.plot(
            data={"sample1": {"x": 1, "y": 2}, "sample2": {"x": 3, "y": 4}},
            headers={"Header1": {"title": "Custom title"}},
            pconfig={
                "name": "Custom table",
                "headers": ["Header1", "Header2"],
                "rows": [["Row1", "Row2"]],
            },
        ),
    )
    report.modules = [module]
    # Should not error:
    multiqc.write_report(force=True, output_dir=str(tmp_path), make_data_dir=False, make_report=False)


def test_software_versions_section(data_dir, tmp_path, capsys):
    multiqc.parse_logs(data_dir / "modules/fastp")
    multiqc.parse_logs(data_dir / "modules/bcftools")
    multiqc.write_report(filename="stdout")  # triggers adding software_versions module
    assert multiqc.list_modules() == ["fastp", "bcftools", "multiqc_software_versions"]


def test_write_report_multiple_times(data_dir, tmp_path):
    multiqc.parse_logs(data_dir / "modules/fastp")
    multiqc.write_report(output_dir=str(tmp_path))
    assert multiqc.list_modules() == ["fastp", "multiqc_software_versions"]

    multiqc.parse_logs(data_dir / "modules/bcftools")
    multiqc.write_report(output_dir=str(tmp_path))
    assert multiqc.list_modules() == ["fastp", "bcftools", "multiqc_software_versions"]


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


def test_user_config(tmp_path, capsys):
    import multiqc

    config_yml = tmp_path / "custom_config.yml"
    expected_title = "Custom Title"
    config_yml.write_text(f'title: "{expected_title}"')

    data_file = tmp_path / "data_mqc.txt"
    data_file.write_text("sample1\t100\nsample2\t200\n")

    multiqc.load_config(config_yml)
    assert multiqc.config.title == expected_title

    multiqc.parse_logs(data_file)
    assert multiqc.config.title == expected_title

    multiqc.write_report(output_dir=str(tmp_path / "output"))
    assert multiqc.config.title == expected_title

    config_yml2 = tmp_path / "custom_config2.yml"
    config_yml2.write_text("""table_cond_formatting_rules:
  column:
    pass:
      - gt: 50
    """)
    multiqc.load_config(config_yml2)
    assert multiqc.config.title == expected_title
    assert multiqc.config.table_cond_formatting_rules["column"]["pass"] == [{"gt": 50}]

    multiqc.parse_logs(data_file)
    assert multiqc.config.title == expected_title
    assert multiqc.config.table_cond_formatting_rules["column"]["pass"] == [{"gt": 50}]

    multiqc.write_report(output_dir=str(tmp_path / "output"))
    assert multiqc.config.title == expected_title
    assert multiqc.config.table_cond_formatting_rules["column"]["pass"] == [{"gt": 50}]
