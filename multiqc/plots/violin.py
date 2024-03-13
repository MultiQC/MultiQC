import logging
from typing import List, Dict, Optional, Union

from multiqc.plots import table_object
from multiqc.utils import config
from multiqc.plots.plotly import violin

logger = logging.getLogger(__name__)


# Load the template so that we can access its configuration
# Do this lazily to mitigate import-spaghetti when running unit tests
_template_mod = None


def get_template_mod():
    global _template_mod
    if not _template_mod:
        _template_mod = config.avail_templates[config.template].load()
    return _template_mod


def plot(data: Union[List[Dict], Dict], headers: Optional[Union[List[Dict], Dict]] = None, pconfig=None):
    """Helper HTML for a violin plot.
    :param data: A list of data dicts
    :param headers: A list of dicts with information
                    for the series, such as colour scales, min and
                    max values etc.
    :param pconfig: plot config dict
    :return: HTML string
    """
    if pconfig is None:
        pconfig = {}
    if not isinstance(data, list):
        data = [data]
    if headers is not None and not isinstance(headers, list):
        headers = [headers]

    # Make datatable objects
    if headers:
        dts = [table_object.DataTable(d, h, pconfig.copy()) for d, h in zip(data, headers)]
    else:
        dts = [table_object.DataTable(d, pconfig=pconfig.copy()) for d in data]

    mod = get_template_mod()
    if "violin" in mod.__dict__ and callable(mod.violin):
        # noinspection PyBroadException
        try:
            return mod.violin(dts)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return violin.plot(dts)
