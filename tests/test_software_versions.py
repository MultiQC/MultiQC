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


def test_software_versions_from_module(data_dir, capsys):
    """
    Verify capturing versions added by another module.
    """
    report.reset()

    mod_dir = data_dir / "modules" / "bismark"
    assert mod_dir.exists() and mod_dir.is_dir()

    multiqc.run(mod_dir, cfg=ClConfig(filename="stdout"))
    captured = capsys.readouterr()
    assert "<td>Bismark</td><td><samp>0.14.0, 0.14.4</samp></td>" in captured.out


def test_software_versions_from_data_and_config(tmp_path, data_dir, capsys):
    """
    Verify finding versions from '*_mqc_versions.yaml' files.
    Verify finding versions from the config section.
    """
    report.reset()

    mod_dir = data_dir / "software_versions"
    assert mod_dir.exists() and mod_dir.is_dir()

    c = {
        "software_versions": {
            "example_module": {
                "example_tool": ["1.4.2"],
            }
        }
    }
    conf_path = tmp_path / "multiqc_config.yaml"
    with open(conf_path, "w") as f:
        yaml.dump(c, f)

    multiqc.parse_logs(mod_dir, config_files=[conf_path])
    multiqc.write_report(filename="stdout")
    captured = capsys.readouterr()
    assert "<td>quast</td><td><samp>4.5.1, 5.1.5</samp></td>" in captured.out
    assert "<td>example_tool</td><td><samp>1.4.2</samp></td>" in captured.out
