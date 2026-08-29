import pytest

from multiqc import report, config, validation
from multiqc.utils import testing


@pytest.fixture
def data_dir():
    return testing.data_dir()


@pytest.fixture(params=["plotly", "echarts"])
def plotting_engine(request):
    """
    Parametrizes a test over both plotting backends by setting `config.plotting_engine`.
    The autouse `reset` fixture below restores `config` to its default after the test.
    """
    config.plotting_engine = request.param
    return request.param


@pytest.fixture(autouse=True)
def reset():
    """
    Reset MultiQC session after use: reset config and report
    """

    yield

    report.reset()
    config.reset()
    validation.reset()
