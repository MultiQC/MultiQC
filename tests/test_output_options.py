import os

import pytest

from multiqc import report, BaseMultiqcModule, write_report


@pytest.fixture()
def stub_modules():
    """
    Set stub modules to make write_report work
    """
    report.modules = [BaseMultiqcModule()]


@pytest.mark.parametrize(
    "options,expected_files",
    [
        ({}, {"multiqc_report.html", "multiqc_data"}),
        ({"filename": "NAME"}, {"NAME.html", "NAME_data"}),
        ({"filename": "NAME.html"}, {"NAME.html", "NAME_data"}),
        ({"title": "My Title"}, {"My-Title_multiqc_report.html", "My-Title_multiqc_report_data"}),
        ({"make_data_dir": False}, {"multiqc_report.html"}),
        ({"make_report": False}, {"multiqc_data"}),
        ({"make_report": False, "make_data_dir": False}, set()),
    ],
)
def test_filename(stub_modules, tmp_path, options, expected_files):
    """
    Verify that the filename option works
    """
    write_report(output_dir=tmp_path, **options)

    assert set(os.listdir(tmp_path)) == expected_files


def test_write_stdout(stub_modules, tmp_path, capsys):
    """
    Verify stdout option: stdout contains only HTML, nothing else is written to disk
    """
    files_before = set(os.listdir(tmp_path))
    os.chdir(tmp_path)

    write_report(filename="stdout")

    captured = capsys.readouterr()
    assert captured.out.startswith("""<!DOCTYPE html>\n<html lang="en">\n<head>""")

    files_after = set(os.listdir(tmp_path))
    assert files_before == files_after
