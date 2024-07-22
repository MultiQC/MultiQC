"""
Generic test for each module: makes sure each file in test-data/data/modules
is found by module and reflected in module.saved_raw_data

Modules can override this test by providing more specific actions for specific files:
e.g. if a file should be skipped, or cause a runtime error, or the raw data should
look in a specific way.
"""

from typing import Callable, List, Union

import pytest

from multiqc import BaseMultiqcModule, config, report, reset, multiqc
from multiqc.core.update_config import update_config

modules = [(k, entry_point) for k, entry_point in config.avail_modules.items() if k != "custom_content"]


@pytest.fixture(scope="module")
def multiqc_reset():
    reset()


@pytest.mark.parametrize("module_id,entry_point", modules)
def test_all_modules(module_id, entry_point, data_dir):
    """
    Verify that all modules do at least something
    """
    config.strict = True

    mod_dir = data_dir / "modules" / module_id
    assert mod_dir.exists() and mod_dir.is_dir()
    report.analysis_files = [mod_dir]
    report.search_files([module_id])

    module_cls: Callable[[], Union[BaseMultiqcModule, List[BaseMultiqcModule]]] = entry_point.load()
    _module = module_cls()
    for m in _module if isinstance(_module, List) else [_module]:
        assert len(report.general_stats_data) > 0 or len(m.sections) > 0


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
        cfg=multiqc.ClConfig(
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
