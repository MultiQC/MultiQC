#!/usr/bin/env python3

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
    filenames = {f[0] for f in report.searchfiles}
    assert filenames == expected


def test_excluded_dirs(init_config, ignored_dirs, search_files):
    """
    Tests that ignored folder names are ignored
    """
    expected_files = {"should_be_included"}
    filenames = {f[0] for f in report.searchfiles}
    assert filenames == expected_files


def test_excluded_paths(init_config, ignored_paths, search_files):
    """
    Tests that ignored *folder* paths are ignored
    """
    expected_files = {"should_be_included"}
    filenames = {f[0] for f in report.searchfiles}
    assert filenames == expected_files
