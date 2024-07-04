"""
Test multiqc running in the command line
"""

import os
import subprocess

import pytest


@pytest.fixture()
def inp_dir(data_dir):
    inp_dir = data_dir / "modules" / "kallisto"
    assert inp_dir.exists() and inp_dir.is_dir()
    return inp_dir


def test_default_commandline(inp_dir, tmp_path):
    """
    Verify that calling `multiqc .` in command line works, returns code 1 and generates the report and data folders
    """
    subprocess.run(["multiqc", inp_dir], cwd=tmp_path, check=True)

    assert set(os.listdir(tmp_path)) == {"multiqc_report.html", "multiqc_data"}


@pytest.mark.parametrize("clean_up", [True, False])
def test_tmpdir_envvar(inp_dir, tmp_path, clean_up):
    tmp_dir = tmp_path / "tmp"
    tmp_dir.mkdir()
    os.environ["TMPDIR"] = str(tmp_dir)

    cwd_dir = tmp_path / "cwd"
    cwd_dir.mkdir()

    cmd = ["multiqc", inp_dir]
    if not clean_up:
        cmd.append("--no-clean-up")
    subprocess.run(cmd, cwd=cwd_dir, check=True)

    assert set(os.listdir(cwd_dir)) == {"multiqc_report.html", "multiqc_data"}
    if clean_up:
        assert set(os.listdir(tmp_dir)) == set()
    else:
        assert set(os.listdir(tmp_dir)) != set()


def test_with_title(inp_dir, tmp_path):
    title = "My Title"
    subprocess.run(["multiqc", inp_dir, "--title", title], cwd=tmp_path, check=True)

    assert set(os.listdir(tmp_path)) == {
        "My-Title_multiqc_report.html",
        "My-Title_multiqc_report_data",
    }


def test_wrong_module(tmp_path):
    """
    Verify that an error is raised when a module is unknown
    """
    result = subprocess.run(
        ["multiqc", tmp_path, "--no-ansi", "-m", "not_a_module"],
        cwd=tmp_path,
        check=False,
        capture_output=True,
    )
    assert result.returncode != 0
    err = result.stderr.decode()
    assert "Invalid value for '-m' / '--module': 'not_a_module' is not one of" in err


def test_unknown_command_line_option(tmp_path):
    result = subprocess.run(
        ["multiqc", tmp_path, "--no-ansi", "--unknown-command-line-flag"],
        cwd=tmp_path,
        check=False,
        capture_output=True,
    )
    assert result.returncode != 0
    err = result.stderr.decode()
    assert "No such option: --unknown-command-line-flag" in err


def test_no_files_found(tmp_path):
    """
    Verify that an error is raised when no files are found
    """
    result = subprocess.run(
        ["multiqc", "--no-ansi", tmp_path],
        cwd=tmp_path,
        check=False,
        capture_output=True,
    )
    assert result.returncode == 0
    err = result.stderr.decode()
    assert "No analysis results found" in err
