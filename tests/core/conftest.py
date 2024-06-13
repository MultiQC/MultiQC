from pathlib import Path

import pytest


@pytest.fixture
def data_dir():
    test_data_dir = Path(__file__).parent.parent.parent / "test-data"
    if not test_data_dir.exists():
        raise FileNotFoundError(
            f"The test data directory expected to be found at {test_data_dir}. Please, "
            f"clone the repository with the test data by changing to the MultiQC repo root, "
            f"and running: `git clone https://github.com/MultiQC/test-data`"
        )

    yield test_data_dir / "data"
