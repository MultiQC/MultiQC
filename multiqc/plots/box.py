"""MultiQC functions to plot a box plot"""

from typing import List, Dict, Union, OrderedDict

import logging

from multiqc import config
from multiqc.plots.plotly.box import BoxT, BoxPlotConfig
from multiqc.plots.plotly import box

logger = logging.getLogger(__name__)

# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod = None


def get_template_mod():
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    return _template_mod


def plot(
    list_of_data_by_sample: Union[Dict[str, BoxT], List[Dict[str, BoxT]]],
    pconfig: Union[Dict, BoxPlotConfig, None],
) -> Union[str, box.BoxPlot]:
    """
    Plot a box plot. Expects either:
    - a dict mapping sample names to data point lists or dicts,
    - a dict mapping sample names to a dict of statistics (e.g. {min, max, median, mean, std, q1, q3 etc.})
    """
    pconf: BoxPlotConfig = BoxPlotConfig.from_pconfig_dict(pconfig)

    # Given one dataset - turn it into a list
    if not isinstance(list_of_data_by_sample, list):
        list_of_data_by_sample = [list_of_data_by_sample]

    for i in range(len(list_of_data_by_sample)):
        if isinstance(list_of_data_by_sample[0], OrderedDict):
            # Legacy: users assumed that passing an OrderedDict indicates that we
            # want to keep the sample order https://github.com/MultiQC/MultiQC/issues/2204
            pass
        elif pconf.sort_samples:
            samples = sorted(list(list_of_data_by_sample[0].keys()))
            list_of_data_by_sample[i] = {s: list_of_data_by_sample[i][s] for s in samples}

    # Make a plot - custom, interactive or flat
    mod = get_template_mod()
    if "box" in mod.__dict__ and callable(mod.bargraph):
        # noinspection PyBroadException
        try:
            return mod.box(list_of_data_by_sample, pconfig)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return box.plot(list_of_data_by_sample, pconf)
