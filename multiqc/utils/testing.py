import os
from pathlib import Path

from multiqc import config


def data_dir():
    test_data_dir = Path(os.environ.get("MULTIQC_TEST_DATA_DIR", config.REPO_DIR / "test-data"))
    if not test_data_dir.exists():
        raise FileNotFoundError(
            f"The test data directory expected to be found at {test_data_dir}. Please, "
            f"clone the repository with the test data by changing to the MultiQC repo root, "
            f"and running: `git clone https://github.com/MultiQC/test-data`"
        )

    return test_data_dir / "data"
