"""
Generic test for each module: makes sure each file in test-data/data/modules
is found by module and reflected in module.saved_raw_data

Modules can override this test by providing more specific actions for specific files:
e.g. if a file should be skipped, or cause a runtime error, or the raw data should
look in a specific way.
"""

from typing import Callable, List, Union

import pytest

from multiqc import BaseMultiqcModule, config, report


modules = [(k, entry_point) for k, entry_point in config.avail_modules.items() if k != "custom_content"]


@pytest.mark.parametrize("module_id,entry_point", modules)
def test_all_modules(module_id, entry_point, data_dir):
    """
    Verify that all modules do at least something
    """
    report.reset()

    mod_dir = data_dir / "modules" / module_id
    assert mod_dir.exists() and mod_dir.is_dir()
    report.analysis_files = [mod_dir]
    report.search_files([module_id])

    module_cls: Callable[[], Union[BaseMultiqcModule, List[BaseMultiqcModule]]] = entry_point.load()
    module = module_cls()

    assert len(report.general_stats_data) > 0 or len(module.sections) > 0
