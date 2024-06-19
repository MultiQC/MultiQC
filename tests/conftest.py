import pytest

from multiqc.utils import testing


@pytest.fixture
def data_dir():
    return testing.data_dir()
