import logging
from typing import Dict, List, Union


from multiqc.templates.plotly.plots.plot import Plot

logger = logging.getLogger(__name__)


# {"name": "SAMPLE1", "color": "#111111", "data": [[x, y], [x, y], ...]}
SampleLineT = Dict[str, Union[str, List]]


def plot(
    datasets: List[List[SampleLineT]],
    pconfig: Dict,
) -> str:
    """
    Build and add the plot data to the report, return an HTML wrapper.
    :param datasets: each dataset is a 2D dict, first keys as sample names, then x:y data pairs
    :param pconfig: dict with config key:value pairs. See CONTRIBUTING.md
    :return: HTML with JS, ready to be inserted into the page
    """
    from multiqc.utils import report

    p = ScatterPlot(pconfig, len(datasets))

    return p.add_to_report(datasets, report)


class ScatterPlot(Plot):
    def __init__(self, pconfig: Dict, *args):
        super().__init__("scatter", pconfig, *args)
