"""Pytest fixtures for the OrthoFinder module tests."""

from pathlib import Path

import pytest


@pytest.fixture
def fixtures_dir():
    """In-repo OrthoFinder output, so the tests need no test-data clone.

    `run_alpha/Comparative_Genomics_Statistics/` holds a three-species run with both
    summary tables, including the trailing distribution blocks that the parsers must
    stop before.
    """
    return Path(__file__).parent / "fixtures"
