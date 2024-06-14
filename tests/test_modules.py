"""
Generic test for each module: makes sure each file in test-data/data/modules
is found by module and reflected in module.saved_raw_data

Modules can override this test by providing more specific actions for specific files:
e.g. if a file should be skipped, or cause a runtime error, or the raw data should
look in a specific way.
"""

from typing import Callable, List, Union

import pytest

import multiqc
from multiqc import BaseMultiqcModule, config, report


@pytest.mark.parametrize("module_id,entry_point", list(config.avail_modules.items())[:1])
def test_all_modules(module_id, entry_point, data_dir):
    mod_dir = data_dir / "modules" / module_id
    assert mod_dir.exists() and mod_dir.is_dir()
    report.analysis_files = [mod_dir]
    report.search_files([module_id])

    module_cls: Callable[[], Union[BaseMultiqcModule, List[BaseMultiqcModule]]] = entry_point.load()
    module = module_cls()
    assert len(module.saved_raw_data) > 0


@pytest.mark.parametrize("module_id,entry_point", list(config.avail_modules.items())[:1])
def test_all_modules_interactive(module_id, entry_point, data_dir):
    mod_dir = data_dir / "modules" / module_id
    assert mod_dir.exists() and mod_dir.is_dir()
    multiqc.parse_logs(mod_dir, run_modules=[module_id])
    assert len(report.saved_raw_data) > 0
    key = list(report.saved_raw_data)[0]
    assert len(report.saved_raw_data[key]) > 0
