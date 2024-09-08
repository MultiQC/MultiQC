"""Test that the sample cleaning logic works as expected."""

import pytest
from multiqc import config
from multiqc.base_module import BaseMultiqcModule


@pytest.fixture
def base_module():
    return BaseMultiqcModule()


def test_no_trim(base_module):
    config.fn_clean_exts[:] = []
    config.fn_clean_trim[:] = []
    assert base_module.clean_s_name("foo.bar.fastq.gz") == "foo.bar.fastq.gz"


def test_default_trim(base_module):
    assert base_module.clean_s_name("foo.bar.fastq.gz") == "foo.bar"


def test_custom_clean_ext(base_module):
    config.fn_clean_exts = ["chop_this_off"]
    assert base_module.clean_s_name("foo.chop_this_off") == "foo"


def test_regex_keep(base_module):
    config.fn_clean_exts = [{"type": "regex_keep", "pattern": "abc..X"}]
    assert base_module.clean_s_name("foo_abc12X_bar") == "abc12X"
    assert base_module.clean_s_name("foo_abc123_bar") == "foo_abc123_bar"
