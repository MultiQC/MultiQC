"""Test that the sample cleaning logic works as expected."""

import pytest
from copy import copy
from multiqc import config
from multiqc.base_module import BaseMultiqcModule


@pytest.fixture
def base_module():
    return BaseMultiqcModule()


@pytest.fixture
def clean_config():
    original_exts = copy(config.fn_clean_exts)
    original_trim = copy(config.fn_clean_trim)

    # This will execute before each test
    yield

    # Restore configuration after each test
    config.fn_clean_exts = original_exts
    config.fn_clean_trim = original_trim


def test_noop(base_module, clean_config):
    # Assume no operation needed
    config.fn_clean_exts[:] = []
    config.fn_clean_trim[:] = []
    assert base_module.clean_s_name("foo.bar.fastq.gz") == "foo.bar.fastq.gz"


def test_default_trim(base_module, clean_config):
    assert base_module.clean_s_name("foo.bar.fastq.gz") == "foo.bar"


def test_custom_clean_ext(base_module, clean_config):
    config.fn_clean_exts = ["chop_this_off"]
    assert base_module.clean_s_name("foo.chop_this_off") == "foo"


def test_regex_keep(base_module, clean_config):
    config.fn_clean_exts = [{"type": "regex_keep", "pattern": "abc..X"}]
    assert base_module.clean_s_name("foo_abc12X_bar") == "abc12X"
    assert base_module.clean_s_name("foo_abc123_bar") == "foo_abc123_bar"
