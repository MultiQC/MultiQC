import os
import tempfile

import pytest

import multiqc
from multiqc import config, report
from multiqc.core.update_config import ClConfig
from multiqc.utils import testing


@pytest.fixture()
def inp_file(tmp_path):
    inp_file = tmp_path / "htseq.tsv"
    inp_file.write_text("""\
feature\tcount
GENE1\t100
GENE2\t200
__no_feature\t413
__ambiguous\t1279
__too_low_aQual\t0
__not_aligned\t3085
__alignment_not_unique\t966""")
    cwd_dir = tmp_path / "work"
    cwd_dir.mkdir()
    os.chdir(cwd_dir)
    return inp_file


def test_write_default(inp_file, tmp_path):
    """
    Verify HTML and data directory with default names are written to current dir
    """
    files_before = set(os.listdir(os.getcwd()))

    multiqc.run(inp_file, cfg=ClConfig(run_modules=["htseq"]))

    files_after = set(os.listdir(os.getcwd()))
    assert files_after - files_before == {"multiqc_report.html", "multiqc_data"}


def test_write_stdout(inp_file, tmp_path, capsys):
    """
    Verify stdout option: stdout contains only HTML, nothing else is written to disk
    """
    files_before = set(os.listdir(tmp_path))
    os.chdir(tmp_path)

    multiqc.run(inp_file, cfg=ClConfig(run_modules=["htseq"], filename="stdout"))

    captured = capsys.readouterr()
    assert captured.out.startswith("""<!DOCTYPE html>\n<html lang="en">\n<head>""")

    files_after = set(os.listdir(tmp_path))
    assert files_before == files_after


@pytest.mark.parametrize("clean_up", [True, False])
def test_no_analysis_found(reset, monkeypatch, tmp_path, clean_up):
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


@pytest.mark.parametrize(
    "options,expected_files",
    [
        ({"filename": "NAME"}, {"NAME.html", "NAME_data"}),
        ({"filename": "NAME.html"}, {"NAME.html", "NAME_data"}),
        ({"title": "My Title"}, {"My-Title_multiqc_report.html", "My-Title_multiqc_report_data"}),
        ({"make_data_dir": False}, {"multiqc_report.html"}),
    ],
)
def test_filename(inp_file, tmp_path, options, expected_files):
    """
    Verify that the filename option works
    """
    multiqc.run(inp_file, cfg=ClConfig(run_modules=["htseq"], **options))

    assert set(os.listdir(os.getcwd())) == expected_files


def test_special_cases(data_dir, tmp_path):
    multiqc.run(
        testing.data_dir() / "special_cases",
        cfg=ClConfig(
            strict=True,
            output_dir=tmp_path,
        ),
    )
    assert len(report.general_stats_data) > 0 or sum(len(m.sections) for m in report.modules) > 0
