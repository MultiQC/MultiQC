"""
Tests for discovering and excluding files
"""

import pytest

from multiqc import config
from multiqc import report
from multiqc.core.file_search import file_search
from multiqc.core.update_config import update_config


@pytest.fixture
def init_config():
    update_config()
    config.sp = {
        "module_included": {
            "fn": "*included.txt",
            "shared": True,
        },
        "module_ignored": {
            "fn": "*ignored.txt",
            "shared": True,
        },
        "module_ignored_gz": {
            "fn": "*ignored.txt.gz",
            "shared": True,
        },
        "module_included_gz": {
            "fn": "*included.txt.gz",
            "shared": True,
        },
        "module_matching_content": {
            "fn": "*.txt",
            "contents": "matching_content",
            "shared": True,
        },
        "module_matching_content_gz": {
            "fn": "*.txt.gz",
            "contents": "matching_content",
            "shared": True,
        },
        "module_filelink": [
            {"fn": "filelink"},
            {"fn": "nested"},
            {"fn": "file"},
        ],
    }
    config.run_modules = list(config.sp.keys())
    config.avail_modules = dict.fromkeys(config.run_modules)


@pytest.fixture
def search_files():
    report.reset_file_search()
    file_search()


@pytest.fixture
def ignored_dirs(data_dir):
    config.analysis_dir = [data_dir / "exclusions/ignore_dirs"]
    config.fn_ignore_dirs = ["ignored", "*_ignored"]


@pytest.fixture
def ignored_paths(data_dir):
    config.analysis_dir = [data_dir / "exclusions/ignore_paths"]
    config.fn_ignore_paths = ["*/*_ignored"]


@pytest.fixture
def ignored_files(data_dir):
    config.analysis_dir = [data_dir / "exclusions/ignore_files"]
    config.fn_ignore_files = ["*.txt.gz"]


@pytest.fixture
def ignore_links(data_dir, request):
    config.analysis_dir = [data_dir / "special_cases/symlinks/linked"]
    config.ignore_symlinks = request.param


@pytest.mark.parametrize(
    ["ignore_links", "expected"],
    [(True, {"file"}), (False, {"filelink", "nested", "file"})],
    indirect=["ignore_links"],
)
def test_symlinked_files_found(init_config, ignore_links, search_files, expected):
    """
    Tests that symlinked files are discovered and ignored properly.
    """
    _assert_expected_files(expected)


def _assert_expected_files(expected_files):
    filenames = set()
    for m, files in report.files.items():
        for f in files:
            filenames.add(f["fn"])
    assert filenames == expected_files


def test_excluded_dirs(init_config, ignored_dirs, search_files):
    """
    Tests that ignored folder names are ignored
    """
    _assert_expected_files({"included.txt"})


def test_excluded_paths(init_config, ignored_paths, search_files):
    """
    Tests that ignored *folder* paths are ignored
    """
    _assert_expected_files({"included.txt"})


def test_fn_ignore_files(init_config, ignored_files, search_files):
    """
    Test that ignored files are still checked against file name if they are in search patterns
    e.g. see Glimpse search pattern that unignores some .txt.gz files
    """
    _assert_expected_files(
        {
            "included.txt",
            "included.txt.gz",
            "matching_content.txt",
        }
    )
