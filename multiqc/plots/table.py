import logging
from typing import List, Dict, Union, Optional

from multiqc.plots import table_object
from multiqc.utils import config
from multiqc.plots.plotly import table

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
    """Return HTML for a MultiQC table.
    :param data: 2D dict, first keys as sample names, then x:y data pairs
    :param headers: list of optional dicts with column config in key:value pairs.
    :param pconfig: plot config dict
    :return: HTML ready to be inserted into the page
    """
    # Make a datatable object
    dt = table_object.DataTable(data, headers, pconfig)

    mod = get_template_mod()
    if "table" in mod.__dict__ and callable(mod.table):
        # Collect unique sample names
        s_names = set()
        for d in dt.data:
            for s_name in d.keys():
                s_names.add(s_name)

        # noinspection PyBroadException
        try:
            return mod.table(dt, s_names, pconfig)
        except:  # noqa: E722
            if config.strict:
                # Crash quickly in the strict mode. This can be helpful for interactive
                # debugging of modules
                raise

    return table.plot([dt])
