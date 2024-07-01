import os

import pytest

import multiqc
from multiqc import report
from multiqc.core.update_config import ClConfig


def test_write_default(data_dir, tmp_path):
    """
    Verify HTML and data directory with default names are written to current dir
    """
    report.reset()
    mod_dir = data_dir / "modules" / "kallisto"
    assert mod_dir.exists() and mod_dir.is_dir()

    files_before = set(os.listdir(tmp_path))
    os.chdir(tmp_path)

    multiqc.run(mod_dir)

    files_after = set(os.listdir(tmp_path))
    assert files_after - files_before == {"multiqc_report.html", "multiqc_data"}


def test_write_stdout(data_dir, tmp_path, capsys):
    """
    Verify stdout option: stdout contains only HTML, nothing else is written to disk
    """
    report.reset()
    mod_dir = data_dir / "modules" / "kallisto"
    assert mod_dir.exists() and mod_dir.is_dir()

    files_before = set(os.listdir(tmp_path))
    os.chdir(tmp_path)

    multiqc.run(mod_dir, cfg=ClConfig(filename="stdout"))

    captured = capsys.readouterr()
    assert captured.out.startswith("""<!DOCTYPE html>\n<html lang="en">\n<head>""")

    files_after = set(os.listdir(tmp_path))
    assert files_before == files_after
