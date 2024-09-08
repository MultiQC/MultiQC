"""
Test the special-case software_versions module.

1. Whether it reports the software version added by another - standard - module via self.add_software_version)(.
2. Whether it from the config section `config.software_versions`
3. Whether it reports versions from found '*_mqc_versions.yaml' files.
"""

import yaml

import multiqc
from multiqc import report
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
                    "example_module": {
                        "example_tool": ["1.4.2"],
                    }
                }
            },
            f,
        )

    # Need some data to be passed besides the bare config
    multiqc.parse_logs(data_dir / "software_versions", config_files=[conf_path])
    multiqc.write_report(filename="stdout")
    captured = capsys.readouterr()
    assert "<td>example_tool</td><td><samp>1.4.2</samp></td>" in captured.out


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
