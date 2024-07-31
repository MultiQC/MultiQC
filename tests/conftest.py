import pytest

from multiqc import report, config
from multiqc.utils import testing


@pytest.fixture
def data_dir():
    return testing.data_dir()


@pytest.fixture(autouse=True)
def reset():
    """
    Reset MultiQC session after use: reset config and report
    """

    yield

    report.reset()
    config.reset()
