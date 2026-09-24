"""Tests for MULTIQC_* environment variable config handling."""

import os

import pytest

from multiqc import config


class TestEnvVarsConfig:
    """Tests for _env_vars_config() parsing MULTIQC_* environment variables."""

    def test_bool_true(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_FORCE", "true")
        config.load_defaults()
        result = config._env_vars_config()
        assert result["force"] is True
        assert isinstance(result["force"], bool)

    def test_bool_false(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_FORCE", "false")
        config.load_defaults()
        result = config._env_vars_config()
        assert result["force"] is False
        assert isinstance(result["force"], bool)

    def test_bool_yes_no(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_FORCE", "yes")
        config.load_defaults()
        result = config._env_vars_config()
        assert result["force"] is True

        monkeypatch.setenv("MULTIQC_FORCE", "no")
        config.load_defaults()
        result = config._env_vars_config()
        assert result["force"] is False

    def test_bool_numeric(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_FORCE", "1")
        config.load_defaults()
        result = config._env_vars_config()
        assert result["force"] is True

        monkeypatch.setenv("MULTIQC_FORCE", "0")
        config.load_defaults()
        result = config._env_vars_config()
        assert result["force"] is False

    def test_make_data_dir_false(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_MAKE_DATA_DIR", "false")
        config.load_defaults()
        result = config._env_vars_config()
        assert result["make_data_dir"] is False
        assert isinstance(result["make_data_dir"], bool)

    def test_make_report_false(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_MAKE_REPORT", "false")
        config.load_defaults()
        result = config._env_vars_config()
        assert result["make_report"] is False
        assert isinstance(result["make_report"], bool)

    def test_int_value(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_PLOTS_FLAT_NUMSERIES", "500")
        config.load_defaults()
        result = config._env_vars_config()
        assert result["plots_flat_numseries"] == 500
        assert isinstance(result["plots_flat_numseries"], int)

    def test_float_value(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_READ_COUNT_MULTIPLIER", "0.001")
        config.load_defaults()
        result = config._env_vars_config()
        assert result["read_count_multiplier"] == 0.001
        assert isinstance(result["read_count_multiplier"], float)

    def test_string_value(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_DATA_FORMAT", "json")
        config.load_defaults()
        result = config._env_vars_config()
        assert result["data_format"] == "json"
        assert isinstance(result["data_format"], str)

    def test_none_typed_value(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_TITLE", "My Report")
        config.load_defaults()
        result = config._env_vars_config()
        assert result["title"] == "My Report"
        assert isinstance(result["title"], str)

    def test_empty_env_var_ignored(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_FORCE", "  ")
        config.load_defaults()
        result = config._env_vars_config()
        assert "force" not in result

    def test_invalid_bool_ignored(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_FORCE", "not_a_bool")
        config.load_defaults()
        result = config._env_vars_config()
        assert "force" not in result

    def test_invalid_int_ignored(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_PLOTS_FLAT_NUMSERIES", "not_a_number")
        config.load_defaults()
        result = config._env_vars_config()
        assert "plots_flat_numseries" not in result

    def test_unknown_key_ignored(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_NONEXISTENT_KEY", "value")
        config.load_defaults()
        result = config._env_vars_config()
        assert "nonexistent_key" not in result

    def test_reserved_name_ignored(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_CONFIG_PATH", "/some/path")
        config.load_defaults()
        result = config._env_vars_config()
        assert "config_path" not in result

    def test_complex_type_ignored(self, monkeypatch):
        monkeypatch.setenv("MULTIQC_FN_IGNORE_DIRS", "some_dir")
        config.load_defaults()
        result = config._env_vars_config()
        assert "fn_ignore_dirs" not in result

    def test_bool_applied_to_config(self, monkeypatch):
        """Verify that parsed boolean env vars are correctly applied to config globals."""
        monkeypatch.setenv("MULTIQC_MAKE_DATA_DIR", "false")
        monkeypatch.setenv("MULTIQC_MAKE_REPORT", "false")
        config.load_defaults()
        env_config = config._env_vars_config()
        config._add_config(env_config)
        assert config.make_data_dir is False
        assert config.make_report is False

    def test_int_applied_to_config(self, monkeypatch):
        """Verify that parsed int env vars are correctly applied to config globals."""
        monkeypatch.setenv("MULTIQC_PLOTS_FLAT_NUMSERIES", "999")
        config.load_defaults()
        env_config = config._env_vars_config()
        config._add_config(env_config)
        assert config.plots_flat_numseries == 999
        assert isinstance(config.plots_flat_numseries, int)
