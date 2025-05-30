import os
from pathlib import Path

import math
import pytest

from multiqc import config, report
from multiqc.modules.samtools import MultiqcModule
from multiqc.utils import testing
from multiqc.utils.testing import BaseModuleTest


@pytest.fixture
def data_dir():
    return testing.data_dir()


@pytest.fixture
def snapshot_config():
    """
    Configuration for snapshot testing.

    Returns a dictionary with common configuration options for snapshot tests.
    """
    return {
        "preserve_module_raw_data": True,
        "strict": True,
        "make_data_dir": False,  # Don't create data directories during testing
    }


def test_data_parsed(data_dir):
    data_subdir = data_dir / "modules/samtools/flagstat"
    for path in os.listdir(data_subdir):
        path = data_subdir / path
        report.analysis_files = [path]
        report.search_files(["samtools"])
        config.preserve_module_raw_data = True
        m = MultiqcModule()
        assert m.saved_raw_data is not None
        assert len(m.saved_raw_data) > 0
        assert m._clean_s_name(Path(path).name) in list(m.saved_raw_data.values())[0]


def slurp_file(data_dir, fname):
    with (data_dir / "modules/samtools/flagstat" / fname).open() as fh:
        return fh.read()


@pytest.fixture
def rep1(data_dir):
    return slurp_file(data_dir, "small.samtools13.flagstat.log.txt")


@pytest.fixture
def rep2(data_dir):
    return slurp_file(data_dir, "small.samtools12.flagstat.log.txt")


def test_rep1(rep1):
    """Test that parsing rep1 produces expected results"""
    from multiqc.modules.samtools.flagstat import parse_single_report

    res1 = parse_single_report(rep1)

    # I expect 13 + 13 + 3 + 3 + 1 things reported in total
    assert len(res1) == 13 + 13 + 3 + 3 + 1

    assert (res1["total_passed"], res1["total_failed"]) == (5414, 0)

    assert res1["flagstat_total"] == 5414

    assert res1["mapped_passed_pct"] == 98.82

    # I expect mapped_failed_pct to be float('nan')
    assert math.isnan(res1["mapped_failed_pct"])


def test_rep2(rep1, rep2):
    """I expect rep2 to give identical results to rep1."""
    from multiqc.modules.samtools.flagstat import parse_single_report

    res1 = parse_single_report(rep1)
    res2 = parse_single_report(rep2)

    # But since nan != nan we have to strip these out.
    nans = [k for k, v in res1.items() if math.isnan(v)]
    for k in nans:
        del res1[k]
        del res2[k]

    assert res1 == res2


# Snapshot testing classes
class TestSamtoolsFlagstatSnapshot(BaseModuleTest):
    """Comprehensive snapshot tests for samtools flagstat module."""

    MODULE_CLASS = MultiqcModule
    MODULE_NAME = "samtools"

    @pytest.fixture
    def flagstat_data_dir(self, data_dir):
        """Get the flagstat-specific data directory."""
        return data_dir / "modules/samtools/flagstat"

    @pytest.fixture
    def module_data_dir(self, data_dir):
        """Get the data directory for this module."""
        return self.get_module_data_dir(data_dir)

    @pytest.fixture
    def module_snapshot(self, data_dir, snapshot_config):
        """Create a snapshot of the module with all available test data."""
        return self.create_module_snapshot(data_dir, snapshot_config)

    @pytest.fixture
    def flagstat_snapshot(self, flagstat_data_dir, snapshot_config):
        """Create a snapshot specifically for flagstat data."""
        # Apply snapshot configuration
        for key, value in snapshot_config.items():
            setattr(config, key, value)

        # Get all flagstat test files
        flagstat_files = list(flagstat_data_dir.glob("*.txt"))

        # Run the module test
        return testing.run_module_test(
            module_class=self.MODULE_CLASS,
            data_files=flagstat_files,
            config_updates=snapshot_config,
        )

    def test_module_data_integrity(self, module_snapshot):
        """Test basic data integrity for the module."""
        self.assert_module_data_integrity(module_snapshot)

    def test_module_complete_snapshot(self, module_snapshot, snapshot):
        """Test the complete module snapshot."""
        complete_snapshot = module_snapshot.get_complete_snapshot()
        assert complete_snapshot == snapshot

    def test_module_raw_data_snapshot(self, module_snapshot, snapshot):
        """Test just the raw data snapshot."""
        raw_data = module_snapshot.get_saved_raw_data()
        assert raw_data == snapshot

    def test_module_general_stats_snapshot(self, module_snapshot, snapshot):
        """Test the general statistics snapshot."""
        general_stats = {
            "data": module_snapshot.get_general_stats_data(),
            "headers": module_snapshot.get_general_stats_headers(),
        }
        assert general_stats == snapshot

    def test_module_sections_snapshot(self, module_snapshot, snapshot):
        """Test the sections snapshot."""
        sections = module_snapshot.get_sections_data()
        assert sections == snapshot

    def test_flagstat_parsing_consistency(self, flagstat_snapshot):
        """Test that flagstat parsing is consistent across different file formats."""
        raw_data = flagstat_snapshot.get_saved_raw_data()

        # Should have parsed data for samtools flagstat
        assert "multiqc_samtools_flagstat" in raw_data
        flagstat_data = raw_data["multiqc_samtools_flagstat"]

        # Check that we have data for the expected samples
        expected_samples = {"small.samtools13.flagstat.log", "small.samtools12.flagstat.log"}
        actual_samples = set(flagstat_data.keys())

        # Should have at least the expected samples (may have more if test data is expanded)
        assert expected_samples.issubset(actual_samples), f"Missing samples: {expected_samples - actual_samples}"

        # All samples should have the same keys (consistent parsing)
        if len(flagstat_data) > 1:
            sample_keys = [set(sample_data.keys()) for sample_data in flagstat_data.values()]
            first_keys = sample_keys[0]
            for i, keys in enumerate(sample_keys[1:], 1):
                assert keys == first_keys, f"Sample {i} has different keys than sample 0"

    def test_flagstat_data_snapshot(self, flagstat_snapshot, snapshot):
        """Snapshot test for flagstat raw data."""
        raw_data = flagstat_snapshot.get_saved_raw_data()
        assert raw_data == snapshot

    def test_flagstat_general_stats_snapshot(self, flagstat_snapshot, snapshot):
        """Snapshot test for flagstat general statistics."""
        general_stats = {
            "data": flagstat_snapshot.get_general_stats_data(),
            "headers": flagstat_snapshot.get_general_stats_headers(),
        }
        assert general_stats == snapshot

    def test_flagstat_sections_snapshot(self, flagstat_snapshot, snapshot):
        """Snapshot test for flagstat sections."""
        sections = flagstat_snapshot.get_sections_data()
        assert sections == snapshot

    def test_flagstat_complete_snapshot(self, flagstat_snapshot, snapshot):
        """Complete snapshot test for flagstat module."""
        complete_snapshot = flagstat_snapshot.get_complete_snapshot()
        assert complete_snapshot == snapshot


# Individual file tests with snapshots
class TestSamtoolsFlagstatIndividualFiles:
    """Test individual flagstat files with snapshots."""

    @pytest.fixture
    def flagstat_data_dir(self, data_dir):
        return data_dir / "modules/samtools/flagstat"

    @pytest.mark.parametrize(
        "filename",
        [
            "small.samtools13.flagstat.log.txt",
            "small.samtools12.flagstat.log.txt",
        ],
    )
    def test_individual_flagstat_file(self, flagstat_data_dir, filename, snapshot, snapshot_config):
        """Test parsing of individual flagstat files."""
        file_path = flagstat_data_dir / filename

        # Apply snapshot configuration
        for key, value in snapshot_config.items():
            setattr(config, key, value)

        # Run the module test on just this file
        module_snapshot = testing.run_module_test(
            module_class=MultiqcModule,
            data_files=[file_path],
            config_updates=snapshot_config,
        )

        # Get the parsed data for this specific file
        raw_data = module_snapshot.get_saved_raw_data()

        # Create a snapshot of just this file's data
        file_snapshot = {
            "filename": filename,
            "parsed_data": raw_data,
            "general_stats": module_snapshot.get_general_stats_data(),
            "sections": module_snapshot.get_sections_data(),
        }

        assert file_snapshot == snapshot
