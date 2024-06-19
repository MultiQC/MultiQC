"""
Tests for discovering and excluding files
"""

import pytest

from multiqc import config
from multiqc import report
from multiqc.core.file_search import file_search
from multiqc.core.update_config import update_config


def _test_search_files(search_patterns, analysis_dir, extra_config, expected_files):
    update_config()

    config.sp = search_patterns
    config.run_modules = list(config.sp.keys())
    config.avail_modules = dict.fromkeys(config.run_modules)

    config.analysis_dir = [analysis_dir]
    if extra_config:
        config.update(extra_config)

    report.reset_file_search()
    file_search()

    filenames = set()
    for m, files in report.files.items():
        for f in files:
            filenames.add(f["fn"])
    assert filenames == expected_files


def test_excluded_dirs(data_dir):
    """
    Tests that ignored folder names are ignored
    """
    _test_search_files(
        search_patterns={
            "module": [
                {"fn": "*included.txt"},
                {"fn": "*ignored.txt"},
            ]
        },
        analysis_dir=data_dir / "exclusions/ignore_dirs",
        extra_config={"fn_ignore_dirs": ["ignored", "*_ignored"]},
        expected_files={"included.txt"},
    )


def test_excluded_paths(data_dir):
    """
    Tests that ignored *folder* paths are ignored
    """
    _test_search_files(
        search_patterns={
            "module": [
                {"fn": "*included.txt"},
                {"fn": "*ignored.txt"},
            ]
        },
        analysis_dir=data_dir / "exclusions/ignore_paths",
        extra_config={"fn_ignore_paths": ["*/*_ignored"]},
        expected_files={"included.txt"},
    )


def test_fn_ignore_files(data_dir):
    """
    Test that ignored files are still checked against file name if they are in search patterns
    e.g. see Glimpse search pattern that unignores some .txt.gz files
    """
    _test_search_files(
        search_patterns={
            "module": [
                {"fn": "*included.txt"},
                {"fn": "*included.txt.gz"},
                {
                    "fn": "*.txt",
                    "contents": "matching_content",
                },
                {
                    "fn": "*.txt.gz",
                    "contents": "matching_content",
                },
            ]
        },
        analysis_dir=data_dir / "exclusions/ignore_files",
        extra_config={"ignored_paths": ["*.txt.gz"]},
        expected_files={
            "included.txt",
            "included.txt.gz",
            "matching_content.txt",
        },
    )


@pytest.mark.parametrize(
    ["ignore_links", "expected_files"],
    [(True, {"file"}), (False, {"filelink", "nested", "file"})],
)
def test_symlinked_files_found(data_dir, ignore_links, expected_files):
    """
    Tests that symlinked files are discovered and ignored properly.
    """
    _test_search_files(
        search_patterns={
            "module_filelink": [
                {"fn": "filelink"},
                {"fn": "nested"},
                {"fn": "file"},
            ],
        },
        analysis_dir=data_dir / "special_cases/symlinks/linked",
        extra_config={"ignore_symlinks": ignore_links},
        expected_files=expected_files,
    )
