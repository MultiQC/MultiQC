import os
import tempfile

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


@pytest.mark.parametrize("clean_up", [True, False])
def test_no_analysis_found(monkeypatch, tmp_path, clean_up):
    """
    Verify that an error is raised when a module is not found
    """
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
