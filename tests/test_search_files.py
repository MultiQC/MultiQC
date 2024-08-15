"""
Tests for discovering and excluding files
"""

from pathlib import Path
from typing import Dict, Set, Union

import pytest
import yaml

from multiqc import config
from multiqc import report
from multiqc.core.exceptions import RunError
from multiqc.core.file_search import file_search


def _test_search_files(
    search_patterns: Dict,
    analysis_dir: Path,
    extra_config: Dict,
    expected_paths_by_module: Dict[str, Set[Union[Path, str]]],
):
    config.sp = search_patterns
    config.run_modules = list(config.sp.keys())
    config.avail_modules = dict.fromkeys(config.run_modules)

    config.analysis_dir = [str(analysis_dir)]
    if extra_config:
        config.update(extra_config)

    file_search()

    expected = {m: set(str(p) for p in paths) for m, paths in expected_paths_by_module.items()}
    found = {m: set(f["fn"] for f in files) for m, files in report.files.items()}
    assert found == expected


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
        analysis_dir=data_dir / "file_search/ignore_dirs",
        extra_config={"fn_ignore_dirs": ["ignored", "*_ignored"]},
        expected_paths_by_module={"module": {"included.txt"}},
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
        analysis_dir=data_dir / "file_search/ignore_paths",
        extra_config={"fn_ignore_paths": ["*/*_ignored"]},
        expected_paths_by_module={"module": {"included.txt"}},
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
        analysis_dir=data_dir / "file_search/ignore_files",
        extra_config={"ignored_paths": ["*.txt.gz"]},
        expected_paths_by_module={
            "module": {
                "included.txt",
                "included.txt.gz",
                "matching_content.txt",
            }
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
            "module": [
                {"fn": "filelink"},
                {"fn": "nested"},
                {"fn": "file"},
            ],
        },
        analysis_dir=data_dir / "special_cases/symlinks/linked",
        extra_config={"ignore_symlinks": ignore_links},
        expected_paths_by_module={"module": expected_files},
    )


def test_filelist(data_dir, tmp_path):
    """
    Test that inputs can be passed and filtered with --file-list
    """
    file1 = tmp_path / "included.txt"
    file2 = tmp_path / "excluded.txt"
    file1.write_text("data1")
    file2.write_text("data2")

    filelist = tmp_path / "filelist.txt"
    filelist.write_text(f"{file1}\n{file2}\nnon_existent.txt\n")

    _test_search_files(
        search_patterns={
            "module": [
                {"fn": "included.txt"},
                {"fn": "non_existent.txt"},
            ]
        },
        analysis_dir=filelist,
        extra_config={"file_list": True},
        expected_paths_by_module={"module": {"included.txt"}},
    )


def test_filelist_all_missing(data_dir, tmp_path):
    """
    Test that if no inputs --file-list found, correct exception is raised
    """
    filelist = tmp_path / "filelist.txt"
    filelist.write_text("non_existent.txt")

    with pytest.raises(RunError):
        _test_search_files(
            search_patterns={
                "module": {"fn": "non_existent.txt"},
            },
            analysis_dir=filelist,
            extra_config={"file_list": True},
            expected_paths_by_module={},
        )


def test_contents_multi_line(tmp_path):
    """
    Test that multi-line contents are correctly matched
    """
    tool1_1 = tmp_path / "tool1-1.txt"
    tool1_1.write_text("""\
metric1: 321
both_tool1_and_tool2_look_for_this_metric: 123
metric3: 321
tool1_also_requires_this_metric: 132
metric5: 321
""")

    tool1_2 = tmp_path / "tool1-2.txt"
    tool1_2.write_text("""\
metric1: 321
tool1_can_alternatively_take_this_metric: 222
""")

    tool2 = tmp_path / "tool2.txt"
    tool2.write_text("""\
metric1: 321
both_tool1_and_tool2_look_for_this_metric: 444
metric3: 321
""")

    sp_patterns: Dict = yaml.safe_load("""
tool1:
  - contents:
    - "both_tool1_and_tool2_look_for_this_metric"
    - "tool1_also_requires_this_metric"
  - contents: "tool1_can_alternatively_take_this_metric"
tool2:
   contents: "both_tool1_and_tool2_look_for_this_metric"
""")

    _test_search_files(
        search_patterns=sp_patterns,
        analysis_dir=tmp_path,
        extra_config={},
        expected_paths_by_module={
            "tool1": {tool1_1.name, tool1_2.name},
            "tool2": {tool2.name},
        },
    )
