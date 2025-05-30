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


# Simplified snapshot testing
class TestSamtoolsFlagstatSnapshot(BaseModuleTest):
    """Snapshot tests for samtools flagstat module."""

    MODULE_CLASS = MultiqcModule
    MODULE_NAME = "samtools"

    def test_flagstat_snapshot(self, data_dir, snapshot, snapshot_config):
        """Test flagstat raw data snapshot."""
        flagstat_data_dir = data_dir / "modules" / "samtools" / "flagstat"
        # Get all flagstat test files
        flagstat_files = list(flagstat_data_dir.glob("*.txt"))
        
        # Run the module test
        module_snapshot = testing.run_module_test(
            module_class=self.MODULE_CLASS,
            data_files=flagstat_files,
            config_updates=snapshot_config,
        )
        
        # Assert basic data integrity
        self.assert_module_data_integrity(module_snapshot)
        
        # Snapshot the raw data output from the module
        raw_data = module_snapshot.get_saved_raw_data()
        assert raw_data == snapshot


# Individual file tests with snapshots
@pytest.mark.parametrize(
    "filename",
    [
        "small.samtools13.flagstat.log.txt",
        "small.samtools12.flagstat.log.txt",
    ],
)
def test_individual_flagstat_file(data_dir, filename, snapshot, snapshot_config):
    """Test parsing of individual flagstat files."""
    flagstat_data_dir = data_dir / "modules" / "samtools" / "flagstat"
    file_path = flagstat_data_dir / filename

    # Run the module test on just this file
    module_snapshot = testing.run_module_test(
        module_class=MultiqcModule,
        data_files=[file_path],
        config_updates=snapshot_config,
    )

    # Just snapshot the raw parsed data for this file
    raw_data = module_snapshot.get_saved_raw_data()
    assert raw_data == snapshot
