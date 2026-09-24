"""Pytest configuration and fixtures for bases2fastq module tests."""

from pathlib import Path

import pytest

from multiqc.utils import testing


@pytest.fixture
def data_dir():
    """Return path to MultiQC test-data repo data directory (test-data/data)."""
    return testing.data_dir()


@pytest.fixture
def fixtures_dir():
    """Return path to in-repo JSON fixtures (no test-data clone required).

    - PairedEndNoProject/RunStats.json (run-level only)
    - PairedEndDefaultProject/RunStats.json, RunManifest.json, Samples/DefaultProject/DefaultProject_RunStats.json
    - PairedEndNoProjectWithLanes/RunStats.json (run-level with Lanes/UnassignedSequences)
    - PairedEndNoProjectLowPolonies/RunStats.json (two samples, one below min_polonies)
    """
    return Path(__file__).parent / "fixtures"
