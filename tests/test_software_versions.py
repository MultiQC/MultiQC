"""
Test the special-case software_versions module.

1. Whether it reports the software version added by another - standard - module via self.add_software_version)(.
2. Whether it from the config section `config.software_versions`
3. Whether it reports versions from found '*_mqc_versions.yaml' files.
"""

import yaml

import multiqc
from multiqc import report
from multiqc.core.special_case_modules.software_versions import MultiqcModule
from multiqc.core.update_config import ClConfig


def test_sorting():
    """
    Versions should be sorted by component
    """
    mod = multiqc.BaseMultiqcModule()

    mod.add_software_version("1.10.1")
    mod.add_software_version("1.2.3")
    mod.add_software_version("1.10")
    mod.add_software_version("1.2.2")

    assert report.software_versions == {
        "base": {
            "base": ["1.2.2", "1.2.3", "1.10", "1.10.1"],
        },
    }


def test_with_software_name():
    """
    Versions should be sorted by component
    """
    mod = multiqc.BaseMultiqcModule()

    mod.add_software_version("1.10.1", software_name="tool1")
    mod.add_software_version("1.2.3", software_name="tool2")
    mod.add_software_version("1.10", software_name="tool2")

    assert report.software_versions == {
        "base": {
            "tool1": ["1.10.1"],
            "tool2": ["1.2.3", "1.10"],
        },
    }


def test_parsing_and_sorting():
    """
    Versions should be parsed for sorting, but represented originally
    """
    mod = multiqc.BaseMultiqcModule()
    versions = [
        "v1.1.1-r505",
        "v2.r505",
        "v3-r505",
        "-r505",
        "r505",
    ]
    for v in versions:
        mod.add_software_version(v)

    assert report.software_versions == {"base": {"base": versions}}


def test_software_versions_from_module(data_dir, capsys):
    """
    Verify versions added by another module are captured
    """
    mod_dir = data_dir / "modules" / "bismark"
    assert mod_dir.exists() and mod_dir.is_dir()

    multiqc.run(mod_dir, cfg=ClConfig(filename="stdout"))
    captured = capsys.readouterr()
    assert "<td>Bismark</td><td><samp>0.14.0, 0.14.4</samp></td>" in captured.out


def test_software_versions_from_config(tmp_path, data_dir, capsys):
    """
    Verify finding versions from the config section.
    """
    conf_path = tmp_path / "multiqc_config.yaml"
    with open(conf_path, "w") as f:
        yaml.dump(
            {
                "software_versions": {
                    "bismark": {
                        "bismark": ["1.4.2"],
                    }
                }
            },
            f,
        )

    # Need some data to be passed besides the bare config
    multiqc.parse_logs(data_dir / "software_versions", config_files=[conf_path])
    multiqc.write_report(filename="stdout")
    captured = capsys.readouterr()
    assert "<td>bismark</td><td><samp>1.4.2</samp></td>" in captured.out


def test_software_versions_from_config_and_module(tmp_path, data_dir, capsys):
    """
    Verify finding versions from the config section and the module.
    """
    mod_dir = data_dir / "modules" / "bismark"
    assert mod_dir.exists() and mod_dir.is_dir()

    multiqc.run(mod_dir, cfg=ClConfig(filename="stdout"))
    captured = capsys.readouterr()

    conf_path = tmp_path / "multiqc_config.yaml"
    with open(conf_path, "w") as f:
        yaml.dump(
            {
                "software_versions": {
                    "bismark": {
                        "bismark": ["1.4.2"],
                    }
                }
            },
            f,
        )

    # Need some data to be passed besides the bare config
    multiqc.parse_logs(data_dir / "software_versions", config_files=[conf_path])
    multiqc.write_report(filename="stdout")
    captured = capsys.readouterr()
    assert "<td>Bismark</td><td><samp>0.14.0, 0.14.4, 1.4.2</samp></td>" in captured.out


def test_software_metadata_defaults_from_module():
    """
    License, license URL and DOI default to the module-level values.
    """
    mod = multiqc.BaseMultiqcModule(
        name="STAR",
        doi="10.1093/bioinformatics/bts635",
        license="MIT License",
        license_url="https://github.com/alexdobin/STAR/blob/master/LICENSE",
    )
    mod.add_software_version("2.7.10a")

    meta = report.software_versions_metadata["STAR"]["STAR"]
    assert meta.license == "MIT License"
    assert meta.license_url == "https://github.com/alexdobin/STAR/blob/master/LICENSE"
    assert meta.doi == ["10.1093/bioinformatics/bts635"]


def test_software_metadata_explicit_per_software():
    """
    Metadata can be provided per software, overriding the module-level values.
    """
    mod = multiqc.BaseMultiqcModule(name="mymodule", license="GPL-3.0")
    mod.add_software_version("1.0.0", software_name="tool1")
    mod.add_software_version("2.0.0", software_name="tool2", license="MIT", doi="10.1000/xyz")

    assert report.software_versions_metadata["mymodule"]["tool1"].license == "GPL-3.0"
    meta2 = report.software_versions_metadata["mymodule"]["tool2"]
    assert meta2.license == "MIT"
    assert meta2.doi == ["10.1000/xyz"]


def test_software_metadata_absent_when_not_provided():
    """
    Modules that don't declare any metadata don't register empty entries.
    """
    mod = multiqc.BaseMultiqcModule(name="plain")
    mod.add_software_version("1.0.0")

    assert "plain" not in report.software_versions_metadata


def test_software_versions_html_license_and_doi_columns():
    """
    The License and DOI columns are rendered only when metadata is present.
    """
    mod = multiqc.BaseMultiqcModule(
        name="STAR",
        doi="10.1093/bioinformatics/bts635",
        license="MIT License",
        license_url="https://github.com/alexdobin/STAR/blob/master/LICENSE",
    )
    mod.add_software_version("2.7.10a")

    html = MultiqcModule._make_versions_html(report.software_versions)
    assert "<th>License</th>" in html
    assert "<th>DOI</th>" in html
    assert '<a href="https://github.com/alexdobin/STAR/blob/master/LICENSE" target="_blank">MIT License</a>' in html
    assert '<a href="https://doi.org/10.1093/bioinformatics/bts635" target="_blank">10.1093/bioinformatics/bts635</a>' in html


def test_software_versions_html_no_extra_columns_without_metadata():
    """
    Without any metadata, the optional columns are not rendered.
    """
    mod = multiqc.BaseMultiqcModule(name="plain")
    mod.add_software_version("1.0.0")

    html = MultiqcModule._make_versions_html(report.software_versions)
    assert "<th>License</th>" not in html
    assert "<th>DOI</th>" not in html


def test_software_versions_from_mqc_files(tmp_path, data_dir, capsys):
    """
    Verify finding versions from '*_mqc_versions.yaml' files.
    """
    mod_dir = data_dir / "software_versions"
    assert mod_dir.exists() and mod_dir.is_dir()

    multiqc.parse_logs(mod_dir)
    multiqc.write_report(filename="stdout")
    captured = capsys.readouterr()
    assert "<td>quast</td><td><samp>4.5.1, 5.1.5</samp></td>" in captured.out
