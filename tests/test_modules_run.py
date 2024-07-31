import tempfile
from pathlib import Path
from typing import Callable, List, Union

import pytest

from multiqc import BaseMultiqcModule, config, report, reset, parse_logs
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.core.update_config import update_config, ClConfig

modules = [(k, entry_point) for k, entry_point in config.avail_modules.items() if k != "custom_content"]


@pytest.fixture(scope="module")
def multiqc_reset():
    reset()


@pytest.mark.parametrize("module_id,entry_point", modules)
def test_all_modules(module_id, entry_point, data_dir):
    """
    Verify that all modules do at least something.

    Modules can add to this test by providing more specific actions for specific files:
    e.g. if a file should be skipped, or cause a runtime error, or the raw data should
    look in a specific way.
    """

    mod_dir = data_dir / "modules" / module_id
    assert mod_dir.exists() and mod_dir.is_dir()

    config.strict = True

    report.analysis_files = [mod_dir]
    report.search_files([module_id])

    module_cls: Callable[[], Union[BaseMultiqcModule, List[BaseMultiqcModule]]] = entry_point.load()
    _module = module_cls()
    for m in _module if isinstance(_module, List) else [_module]:
        assert len(report.general_stats_data) > 0 or len(m.sections) > 0


@pytest.mark.parametrize("module_id,entry_point", modules)
def test_ignore_samples(module_id, entry_point, data_dir):
    """
    Verify all modules call self.ignore_samples() correctly.
    """
    mod_dir = data_dir / "modules" / module_id
    assert mod_dir.exists() and mod_dir.is_dir()

    config.sample_names_ignore = ["*"]

    report.analysis_files = [mod_dir]
    report.search_files([module_id])

    module_cls: Callable[[], Union[BaseMultiqcModule, List[BaseMultiqcModule]]] = entry_point.load()
    with pytest.raises(ModuleNoSamplesFound):
        _module = module_cls()


@pytest.mark.parametrize(
    ["config_options", "expected_to_write"],
    [
        (dict(), True),
        ({"make_data_dir": False}, False),
        ({"filename": "stdout"}, False),
    ],
)
def test_write_data_file(monkeypatch, tmp_path, config_options, expected_to_write):
    """
    Test module.write_data_file() write something
    """
    (tmp_path / "multiqc_tmp").mkdir()
    monkeypatch.setattr(tempfile, "mkdtemp", lambda: tmp_path / "multiqc_tmp")
    config.update(config_options)
    module = BaseMultiqcModule()
    module.write_data_file({"Sample": {"key": "value"}}, "multiqc_mymodule")

    expected_path = tmp_path / "multiqc_tmp" / "multiqc_data" / "multiqc_mymodule.txt"
    if expected_to_write:
        assert expected_path.exists()
        assert expected_path.open().read().strip() == """Sample\tkey\nSample\tvalue""".strip()
    else:
        assert not expected_path.exists()


@pytest.mark.parametrize(
    "use_filename_as_sample_name,fn_clean_sample_names,prepend_dirs,expected_sample_name",
    [
        (None, None, None, "SAMPLE_FROM_CONTENTS"),
        (True, None, None, "SAMPLE_FROM_FILENAME"),
        (None, False, None, "SAMPLE_FROM_CONTENTS.fastq.gz"),
        (True, False, None, "SAMPLE_FROM_FILENAME.stderr"),
        (None, None, True, "subdir | SAMPLE_FROM_CONTENTS"),
        (True, None, True, "subdir | SAMPLE_FROM_FILENAME"),
    ],
)
def test_use_filename_as_sample_name(
    multiqc_reset, tmp_path, use_filename_as_sample_name, fn_clean_sample_names, prepend_dirs, expected_sample_name
):
    """
    Verify that `--fn_as_s_name`, `--fullnames`, and `--dirs` works
    """
    report.reset()

    MODULE_NAME = "trimmomatic"
    (tmp_path / "subdir").mkdir()
    input_file = tmp_path / "subdir" / "SAMPLE_FROM_FILENAME.stderr"
    input_file.write_text("""\
TrimmomaticSE: Started with arguments: SAMPLE_FROM_CONTENTS.fastq.gz
Input Reads: 39733090 Surviving: 32590558 (82.02%) Dropped: 7142532 (17.98%)
TrimmomaticSE: Completed successfully""")

    update_config(
        cfg=ClConfig(
            run_modules=[MODULE_NAME],
            use_filename_as_sample_name=use_filename_as_sample_name,
            fn_clean_sample_names=fn_clean_sample_names,
            prepend_dirs=prepend_dirs,
            dirs_depth=1 if prepend_dirs else None,
        )
    )

    report.analysis_files = [input_file]
    report.search_files([MODULE_NAME])

    from multiqc.modules.trimmomatic.trimmomatic import MultiqcModule

    m = MultiqcModule()

    assert expected_sample_name in m.saved_raw_data[f"multiqc_{MODULE_NAME}"]


def test_path_filters(multiqc_reset, tmp_path, data_dir):
    search_path = data_dir / "modules" / "adapterremoval"
    assert search_path.exists() and search_path.is_dir()

    expected_pe_files = {
        "paired_end_collapsed/pec1.settings",
        "paired_end_collapsed/pec2.settings",
        "paired_end_noncollapsed/penc1.settings",
        "paired_end_noncollapsed/penc2.settings",
    }
    expected_se_files = {
        "single_end/se.settings",
    }

    assert all((search_path / fn).exists() for fn in expected_pe_files)
    assert all((search_path / fn).exists() for fn in expected_pe_files)

    parse_logs(
        search_path,
        module_order=[
            {
                "adapterremoval": {
                    "name": "adapterremoval (single end)",
                    "anchor": "my_anchor_se",
                    "path_filters": ["*/se.*"],
                },
            },
            {
                "adapterremoval": {
                    "name": "adapterremoval (paired end)",
                    "anchor": "my_anchor_pe",
                    "path_filters": ["*/pec?.*", "*/penc?.*"],
                },
            },
        ],
    )

    assert len(report.modules) == 2
    assert len(report.general_stats_data) == 2
    assert report.modules[0].name == "adapterremoval (single end)"
    assert report.modules[1].name == "adapterremoval (paired end)"
    assert report.modules[0].saved_raw_data["multiqc_adapter_removal_my_anchor_se"].keys() == {
        Path(fn).name for fn in expected_se_files
    }
    assert report.modules[1].saved_raw_data["multiqc_adapter_removal_my_anchor_pe"].keys() == {
        Path(fn).name for fn in expected_pe_files
    }
    assert report.general_stats_data[0].keys() == {Path(fn).name for fn in expected_se_files}
    assert report.general_stats_data[1].keys() == {Path(fn).name for fn in expected_pe_files}
