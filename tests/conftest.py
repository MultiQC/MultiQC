import pytest

import multiqc
from multiqc.utils import testing


@pytest.fixture
def data_dir():
    return testing.data_dir()


@pytest.fixture
def reset():
    multiqc.reset()
